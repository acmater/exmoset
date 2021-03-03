import unittest
import numpy as np
import pandas as pd

from molsim.atom import ContainsAtom
from molsim.bond import ContainsBond
from molsim.substructure import Substructure
from molsim.data import *
from molsim.molecule import MolProp

test_mols = molecules4

class TestSubstructure(unittest.TestCase):
    def test_substructures(self):
        substructures = Substructure(test_mols,substructures=["[NH2]"])
        assert sum(substructures["[NH2]"]) == 1, "Incorrectly identifying substructures"


class TestDatabaseExtraction(unittest.TestCase):
    def test_molprop(self):
        df     = pd.read_csv("molsim/data/QM9_Data.csv",index_col="SMILES")
        sub_df = df.loc[test_mols]
        molprops = MolProp(test_mols,properties=["Dipole Moment","ZPVE"],df=sub_df)
        assert len(molprops["Dipole Moment"] == 86), "Molecular Property dictionary is not being properly generated."
        assert max(molprops["ZPVE"] == 0.216981), "Issue in the ZPVE section of the dictionary."

class TestContainsAtom(unittest.TestCase):
    global atom_props
    atom_props = ContainsAtom(test_mols,["C","N","O"])
    def test_contains_C(self):
        assert atom_props.entropy(atom_props["C"]) == 0, "Entropy calculation for contains C is broken"
    def test_contains_O(self):
        assert atom_props.entropy(atom_props["O"])


class TestContainsBond(unittest.TestCase):
    global bond_props
    bond_props = ContainsAtom(test_mols,["SINGLE","DOUBLE","TRIPLE"])
    def test_contains_SINGLE(self):
        assert bond_props.entropy(bond_props["SINGLE"]) == 0, "Entropy calculation for contains Single bond is broken"
    def test_contains_TRIPLE(self):
        assert bond_props.entropy(bond_props["TRIPLE"]) == 0, "Entropy Calculation for contains Triple bond is broken"

if __name__ == "__main__":
    unittest.main()
