"""Check exact positional labels, immutable templates and compiled-copy transfers."""
from collections import defaultdict
import random
import unittest
from rdkit import Chem
from carbon_flow_tools.atom_helpers import labeled_smiles, _labeled_smiles, _structure_template
from carbon_flow_tools.carbon_trace import map_one_reaction, reaction_transfer
from scripts.test_carbon_tracing import reaction


def reference_smiles(smiles, assignment):
    labels={(int(p),int(a)):int(str(label)[1:]) for p,a,label in assignment}
    output=[]
    for part,fragment in enumerate(smiles.split('.')):
        molecule=Chem.MolFromSmiles(fragment)
        if molecule is None:raise ValueError('Invalid SMILES: '+fragment)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(labels.get((part,atom.GetIdx()),0))
        output.append(Chem.MolToSmiles(molecule,canonical=True,isomericSmiles=True))
    return '.'.join(output)


class EnumerationOptimizationTests(unittest.TestCase):
    def test_templates_match_original_writer_without_mutating_shared_structure(self):
        randomizer=random.Random(7301)
        structures=['CCO','C[C@H](O)F','C[C@@H](O)F','C1CCCCC1','c1cc[nH]c1',
            '[13CH3]C(=O)[O-]','[CH3:99][CH2:12]O','*SC(=O)CCC(=O)OC',
            'C.CO','[Na+].CC(=O)[O-]','P(=O)([O-])([O-])OCC',
            '*C(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(*)=O']
        for smiles in structures:
            positions=[(part,atom.GetIdx()) for part,fragment in enumerate(smiles.split('.'))
                       for atom in Chem.MolFromSmiles(fragment).GetAtoms()]
            for _ in range(20):
                assignment=tuple((p,a,'C'+str(randomizer.randrange(1,12))) for p,a in positions if randomizer.random()<.5)
                self.assertEqual(labeled_smiles(smiles,assignment),reference_smiles(smiles,assignment))
            for fragment in smiles.split('.'):
                self.assertTrue(all(a.GetAtomMapNum()==0 for a in _structure_template(fragment).GetAtoms()))
        self.assertEqual(labeled_smiles('CC',((0,-1,'C1'),(5,0,'C2'),(0,99,'C3'))),'CC')
        self.assertEqual(_structure_template.cache_info().maxsize,512)
        self.assertEqual(_labeled_smiles.cache_info().maxsize,8192)

    def test_copy_order_and_required_source_filter_remain_exact(self):
        rxn,components=reaction('2 A[c]','B[c]',{'A[c]':'C','B[c]':'CC'},'[CH4:1].[CH4:2]>>[CH3:1][CH3:2]')
        states={'A[c]':{1:((0,0,'C1'),),2:((0,0,'C2'),)}}
        catalog={}
        produced,rows,next_id=map_one_reaction(rxn,components,states,0,catalog,1,{'B[c]':'CC'},
            {('A[c]',1):0,('A[c]',2):4},required_labels={'C2'})
        self.assertEqual([r['SourceIDs'] for r in rows],['1&2','2&1','2&2'])
        self.assertEqual([r['Step'] for r in rows],[5,5,5])
        self.assertEqual([r['ProductionID'] for r in rows],[1,1,2])
        self.assertEqual(next_id,3)
        self.assertEqual([r['CarbonLabels'] for r in rows],['C1,C2','C1,C2','C2,C2'])
        original=reaction_transfer(rxn,components)
        self.assertIs(reaction_transfer(rxn,components),original)
        rxn.source_text='A[c]'
        # Valid changes to the active stoichiometry invalidate the fixed plan.
        changed=(components[0][:1],components[1])
        self.assertIsNot(reaction_transfer(rxn,changed),original)

    def test_absent_substrate_and_earliest_steps_are_not_approximated(self):
        rxn,components=reaction('A[c] + B[c]','P[c]',{'A[c]':'C','B[c]':'C','P[c]':'CC'},'[CH4:1].[CH4:2]>>[CH3:1][CH3:2]')
        states={'A[c]':{0:((0,0,'C1'),)}}
        first={('A[c]',0):2}
        values=map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CC'},first)[1]
        self.assertEqual(values[0]['SourceIDs'],'0 + 0')
        self.assertEqual(values[0]['CarbonLabels'],'C1')
        states['B[c]']={9:((0,0,'C2'),)};first[('B[c]',9)]=7
        values=map_one_reaction(rxn,components,states,0,{},1,{'P[c]':'CC'},first)[1]
        self.assertEqual(values[0]['SourceIDs'],'0 + 9')
        self.assertEqual(values[0]['Step'],8)
        self.assertEqual(values[0]['CarbonLabels'],'C1,C2')


if __name__=='__main__':unittest.main()
