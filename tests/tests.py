import unittest
import numpy as np
import pandas as pd
from rdkit import Chem

from exmoset.data import *
from exmoset.molecule import Molecule
from exmoset.fingerprints import *
from exmoset.labels import *

test_mols = molecules4

# Is there a better way to do the global portions of the code below? Intuitively I should be able
# to customize the initialization of each class and provide it as an attribute, but that doesn't seem to work.

class TestMolecule(unittest.TestCase):
    def test_mol_gen(self):
        test = Molecule("CCC")
        assert test, "Molecule object not correctly generated"
    def test_mol_gen_rdkit(self):
        mol_converters = {"mol" : Chem.MolFromSmiles("CCC")}
        test2 = Molecule("CCC",**mol_converters)
        assert isinstance(test2.mol,Chem.rdchem.Mol), "Extra attribute generation not working for Molecule class"
    def test_mol_eq(self):
        assert Molecule("CCC",mol=Chem.MolFromSmiles("CCC")) == Molecule("CCC"), "Molecule equality not working"

class TestFingerprint(unittest.TestCase):
    def test_fingerprint(self):
        assert Fingerprint(name="Contains C",
                    context="molecule",
                    label_type="binary",
                    calculator="add",
                    mol_format="smiles"), "Fingerprint generation did not work properly"

class TestBinaryLabel(unittest.TestCase):
    global binary
    binary = Binary("binary",np.ones((100,1)))
    def test_av(self):
        assert binary.av == 1, "Binary label averaging not working."
    def test_entropy(self):
        assert binary.entropy == 0, "Entropy testing for Binary labels not working."
    def test_name(self):
        assert binary.name == "binary", "Binary name not working properly."

class TestMulticlassLabel(unittest.TestCase):
    global multi
    multi = Multiclass("multi",np.arange(10))
    def test_av(self):
        assert multi.av == 4, "Multiclass label averaging not working."
    def test_entropy(self):
        assert multi.entropy != 0, "Entropy testing for Multiclass labels not working."
    def test_name(self):
        assert multi.name == "multi", "Multiclass label name not working properly."

class TestContinuousclassLabel(unittest.TestCase):
    global cont
    cont = Continuous("cont",np.arange(10))
    def test_av(self):
        assert cont.av == 4.5, "Continuous label averaging not working."
    def test_entropy(self):
        assert cont.entropy != 0, "Entropy testing for Continuous labels not working."
    def test_name(self):
        assert cont.name == "cont", "Continuous label name not working properly."

class TestMolSet(unittest.TestCase):
    global analysis
    fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints
    try:
        analysis = MolSet(molecules4,
                        fingerprints = fingerprints,
                        mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                        significance=0.1,
                        file="data/QM9_Data.csv")
    except Exception:
        "Molset Class could not be created"




if __name__ == "__main__":
    unittest.main()
