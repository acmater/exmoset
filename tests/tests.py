import unittest
import numpy as np
import pandas as pd

from molsim.atom import ContainsAtom
from molsim.bond import ContainsBond
from molsim.substructure import Substructure
from molsim.data import *
from molsim.molecule import MolProp

test_mols = molecules4

# Is there a better way to do the global portions of the code below? Intuitively I should be able
# to customize the initialization of each class and provide it as an attribute, but that doesn't seem to work.

class TestFailedMoleculeConversion(unittest.TestCase):
    def test_failed_import(self):
        assert len(Substructure(broken,substructures=["[CC]"])["[CC]"]) == 0, "The system is not removing failed smile conversions"

class TestSubstructure(unittest.TestCase):
    global substructures
    substructures = Substructure(test_mols,substructures=["[NH2]"])
    def test_substructures(self):
        assert sum(substructures["[NH2]"]) == 1, "Incorrectly identifying substructures"

class TestDatabaseExtraction(unittest.TestCase):
    global df, sub_df, molprops
    df     = pd.read_csv("molsim/data/QM9_Data.csv",index_col="SMILES")
    sub_df = df.loc[test_mols]
    molprops = MolProp(test_mols,properties=["Dipole Moment","ZPVE"],df=sub_df)
    def test_molprop(self):
        assert len(molprops["Dipole Moment"]) == 86, "Molecular Property dictionary is not being properly generated."

    def test_max_ZPVE(self):
        assert max(molprops["ZPVE"]) == 0.216981, "Issue in the ZPVE section of the dictionary."

    def test_fake_prop(self):
        with self.assertRaises(KeyError):
            molprops["Electronic Spatial Extent"]

class TestContainsAtom(unittest.TestCase):
    global atom_props
    atom_props = ContainsAtom(test_mols,["C","N","O"])

    def test_contains_C(self):
        assert atom_props.entropy(atom_props["C"]) == 0, "Entropy calculation for contains C is broken"
    def test_contains_O(self):
        assert atom_props.entropy(atom_props["O"]) > 0, "Entropy Calculation for contains O is broken"
    def test_contains_F(self):
        with self.assertRaises(KeyError):
            atom_props["F"]

class TestContainsBond(unittest.TestCase):
    global bond_props
    bond_props = ContainsAtom(test_mols,["SINGLE","DOUBLE","TRIPLE"])

    def test_contains_SINGLE(self):
        assert bond_props.entropy(bond_props["SINGLE"]) == 0, "Entropy calculation for contains Single bond is broken"
    def test_contains_TRIPLE(self):
        assert bond_props.entropy(bond_props["TRIPLE"]) == 0, "Entropy Calculation for contains Triple bond is broken"

if __name__ == "__main__":
    unittest.main()
