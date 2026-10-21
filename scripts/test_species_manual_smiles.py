"""Species/name manual structures are shared by conditions, not unrelated species."""
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from carbon_flow_tools.smiles_completion import complete_smiles, load_manual, species_name

class SpeciesManualTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.root=Path(self.directory.name)
        self.manual=self.root/'manual.csv'
        self.options=dict(manual_file=self.manual,cache_file=self.root/'cache.json',lookup=False,timeout=1,min_interval=0)
    def tearDown(self):
        self.directory.cleanup()
    def write(self,rows):
        with self.manual.open('w',newline='') as handle:
            writer=csv.writer(handle)
            writer.writerow(['Species','MetaboliteName','SMILES'])
            writer.writerows(rows)
    def test_same_species_conditions_compartments_and_different_ids(self):
        self.write([['Escherichia_coli','ethanol','CCO']])
        table=pd.DataFrame([['different_id','ethanol[c]'],['another_id','ethanol[e]']],columns=['MetaboliteID','MetaboliteName'])
        with patch('carbon_flow_tools.smiles_completion.fetch_properties',side_effect=AssertionError('No lookup')):
            for name in ('Escherichia_coli_non_IUPD','Escherichia_coli_IUPD','Escherichia_coli_glucose_xylose'):
                result=complete_smiles(table,name,self.options)
                self.assertTrue(result.SMILES_Manual.eq('CCO').all())
            result=complete_smiles(table,'Saccharomyces_cerevisiae_IUPD',self.options)
            self.assertTrue(result.SMILES_Manual.isna().all())
    def test_source_directory_can_identify_species_for_custom_case_name(self):
        self.assertEqual(species_name('GlucoseOnly',dict(path='data/Escherichia_coli/iML1515.mat')),'Escherichia_coli')
        with self.assertRaises(ValueError):
            species_name('Saccharomyces_cerevisiae_IUPD',dict(path='data/Escherichia_coli/iML1515.mat'))
    def test_conflicting_names_are_rejected_before_choosing_a_structure(self):
        self.write([['Escherichia_coli','compound[c]','CCO'],['Escherichia_coli','compound[e]','COC']])
        with self.assertRaisesRegex(ValueError,'Conflicting'):
            load_manual(self.manual)
    def test_formula_and_charge_conflicts_are_not_overridden(self):
        self.write([['Escherichia_coli','ethanol','CCO']])
        for formula,charge in [('CH4O',0),('C2H6O',-1)]:
            table=pd.DataFrame([dict(MetaboliteID='m1',MetaboliteName='ethanol[c]',MolecularFormula=formula,Charge=charge)])
            with self.assertRaisesRegex(ValueError,'conflict'):
                complete_smiles(table,'Escherichia_coli_non_IUPD',self.options)
    def test_ai_cache_fills_missing_only_and_keeps_review_provenance(self):
        self.write([])
        ai=self.root/'ai_smiles.csv'
        ai.write_text('Species,MetaboliteName,SMILES,Source\nEscherichia_coli,ethanol,CCO,https://example.org/evidence\n')
        table=pd.DataFrame([dict(MetaboliteID='a',MetaboliteName='ethanol[c]',MolecularFormula='C2H6O',Charge=0)])
        result=complete_smiles(table,'Escherichia_coli_IUPD',self.options)
        self.assertEqual(result.iloc[0].SMILES_Retrieved,'CCO')
        self.assertTrue(result.iloc[0].StructureReviewRequired)
        self.assertEqual(result.iloc[0].StructureReference,'https://example.org/evidence')
        self.assertTrue(result.SMILES_Manual.isna().all())
        table['SMILES_Model']='CCO'
        existing=complete_smiles(table,'Escherichia_coli_IUPD',self.options)
        self.assertTrue(existing.SMILES_Retrieved.isna().all())
        self.write([['Escherichia_coli','ethanol','CCO']])
        manual=complete_smiles(table,'Escherichia_coli_IUPD',self.options)
        self.assertEqual(manual.iloc[0].StructureSource,'manual table')
        self.assertFalse(manual.iloc[0].StructureReviewRequired)

    def test_schema_has_exactly_three_fields(self):
        self.manual.write_text('Model,MetaboliteID,MetaboliteName,SMILES\n')
        with self.assertRaisesRegex(ValueError,'exactly'):
            load_manual(self.manual)

if __name__=='__main__':
    unittest.main()