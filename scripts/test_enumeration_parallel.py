"""Verify pure parallel chunks cannot change order, IDs, state or safety stops."""
from copy import deepcopy
import os
import time
import unittest
from unittest.mock import patch

from carbon_flow_tools.enumeration import EnumerationPool, _candidate_batch
from carbon_flow_tools.carbon_trace import map_one_reaction
from carbon_flow_tools.config import MemoryPressureError
from carbon_flow_tools.excel_io import WorkbookSizeLimitError
from carbon_flow_tools.task_scheduler import run_tasks
from scripts.test_carbon_tracing import reaction


def delayed_batch(*args):
    if args[2][0][0][0] == 1:
        time.sleep(.08)
    return _candidate_batch(*args)


def bad_batch(*args):
    raise ValueError('candidate fixture failure')


def dead_batch(*args):
    os._exit(9)


def fixture():
    rxn,components=reaction('2 A[c]','P[c]',{'A[c]':'C','P[c]':'CCO'},
        '[CH4:1].[CH4:2]>>[CH3:1][CH2:2]O')
    states={'A[c]':{i:((0,0,'C'+str(i)),) for i in range(1,9)}}
    first={('A[c]',i):i%3 for i in range(1,9)}
    return rxn,components,states,first


def nested_probe():
    from multiprocessing import current_process
    rxn,components,states,first=fixture()
    pool=EnumerationPool(2,3,2)
    try:
        with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS',0):
            result=map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CCO'},first,candidate_pool=pool)
        return current_process().daemon,result
    finally:
        pool.close()


class ParallelEnumerationTests(unittest.TestCase):
    def test_out_of_order_completion_retains_every_row_and_id(self):
        rxn,components,states,first=fixture()
        before=deepcopy((states,first));serial_catalog={};parallel_catalog={}
        args=(rxn,components,states,0)
        expected=map_one_reaction(*args,catalog=serial_catalog,next_id=1,
            smiles_by_name={'P[c]':'CCO'},first_steps=first,required_labels={'C3'},unordered_sources=True)
        environment=os.environ.get('OPENBLAS_NUM_THREADS')
        pool=EnumerationPool(2,3,2)
        try:
            with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS',0), \
                 patch('carbon_flow_tools.enumeration.CHUNK_SIZE',3), \
                 patch('carbon_flow_tools.enumeration._candidate_batch',delayed_batch):
                actual=map_one_reaction(*args,catalog=parallel_catalog,next_id=1,
                    smiles_by_name={'P[c]':'CCO'},first_steps=first,required_labels={'C3'},
                    unordered_sources=True,candidate_pool=pool)
            self.assertIsNotNone(pool.pool)
            self.assertEqual(actual,expected)
            self.assertEqual(parallel_catalog,serial_catalog)
            self.assertEqual((states,first),before)
        finally:
            pool.close()
        self.assertEqual(os.environ.get('OPENBLAS_NUM_THREADS'),environment)

    def test_prefetched_combinations_are_not_reported_as_completed(self):
        rxn, components, states, first = fixture()
        states['A[c]'] = {i: ((0, 0, 'C' + str(i)),) for i in range(1, 25)}
        seen, counts = [], []
        class TrackingProgress:
            def begin(self, label, total):
                self.total = total
            def update(self, completed):
                # Each combination in this fixture produces exactly one row.
                self_test.assertLessEqual(completed, len(seen))
                counts.append(completed)
        self_test = self
        progress = TrackingProgress()
        pool = EnumerationPool(2, 3, 2)
        try:
            with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS', 0), \
                 patch('carbon_flow_tools.enumeration.CHUNK_SIZE', 128):
                map_one_reaction(rxn, components, states, 1, {}, 1, {'P[c]': 'CCO'},
                    record_sink=seen.append, candidate_pool=pool, progress=progress)
        finally:
            pool.close()
        self.assertEqual(counts, sorted(counts))
        self.assertEqual(counts[-1], 576)
        self.assertEqual(progress.total, 576)
        self.assertEqual(len(seen), 576)

    def test_file_stop_cancels_prefetched_chunks_without_writing_more_rows(self):
        rxn,components,states,first=fixture();seen=[]
        def sink(row):
            seen.append(row)
            if len(seen)==10:raise WorkbookSizeLimitError('fixture size limit')
        pool=EnumerationPool(2,3,2)
        with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS',0), \
             patch('carbon_flow_tools.enumeration.CHUNK_SIZE',3):
            with self.assertRaisesRegex(WorkbookSizeLimitError,'fixture size limit'):
                map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CCO'},first,
                                 record_sink=sink,candidate_pool=pool)
        self.assertEqual(len(seen),10)
        self.assertIsNone(pool.pool)

    def test_worker_failure_remains_explicit_and_closes_pool(self):
        rxn,components,states,first=fixture();pool=EnumerationPool(2,3,2)
        with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS',0), \
             patch('carbon_flow_tools.enumeration._candidate_batch',bad_batch):
            with self.assertRaisesRegex(ValueError,'candidate fixture failure'):
                map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CCO'},first,candidate_pool=pool)
        self.assertIsNone(pool.pool)

    def test_worker_death_is_reported_instead_of_waiting_for_a_lost_task(self):
        rxn,components,states,first=fixture();pool=EnumerationPool(2,3,2)
        with patch('carbon_flow_tools.enumeration.MIN_PARALLEL_COMBINATIONS',0), \
             patch('carbon_flow_tools.enumeration._candidate_batch',dead_batch):
            with self.assertRaisesRegex(RuntimeError,'worker exited unexpectedly'):
                map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CCO'},first,candidate_pool=pool)
        self.assertIsNone(pool.pool)

    def test_outer_analysis_processes_allow_inner_pools(self):
        with patch('carbon_flow_tools.task_scheduler.available_memory',return_value=(32,32)):
            values=run_tasks(nested_probe,[(),()],2,3,allow_children=True,cpu_threads=1)
        self.assertTrue(all(not daemon for daemon,_ in values))
        self.assertEqual(values[0][1],values[1][1])

    def test_invalid_process_setting_does_not_open_output_storage(self):
        from openpyxl import Workbook
        from carbon_flow_tools.pipeline import trace_analysis
        from scripts.test_output_limits import flow_fixture
        import tempfile
        from pathlib import Path
        rxn,metabolites,_=flow_fixture()
        with tempfile.TemporaryDirectory() as folder, \
             patch('carbon_flow_tools.pipeline.FlowOutput') as store:
            with self.assertRaises(ValueError):
                trace_analysis('Invalid',dict(carbon_sources=['A[c]'],max_enumeration_processes=0),
                    [rxn],metabolites,2,Workbook(),Path(folder)/'output.xlsx')
            store.assert_not_called()

    def test_capacity_and_invalid_options(self):
        for invalid in (0,-1,True,1.5):
            with self.assertRaises(ValueError):EnumerationPool(4,4,invalid)
        self.assertEqual(EnumerationPool(1,4).workers,1)
        with patch('carbon_flow_tools.enumeration.available_memory',return_value=(128,128)), \
             patch('carbon_flow_tools.enumeration.available_cpus',return_value=128):
            self.assertEqual(EnumerationPool(32,20).workers,32)
            self.assertEqual(EnumerationPool(128,100).workers,128)
            self.assertEqual(EnumerationPool(32,2).workers,2)
        with patch('carbon_flow_tools.enumeration.available_memory',side_effect=OSError):
            self.assertEqual(EnumerationPool(32,20).workers,1)


if __name__=='__main__':unittest.main()
