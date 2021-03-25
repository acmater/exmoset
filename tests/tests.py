import unittest
import numpy as np
import pandas as pd
from rdkit import Chem

from exmoset.data import *
from exmoset.molecule import Molecule
from fingerprints import *
from exmoset import MolSpace

import matplotlib
import matplotlib.pyplot as plt

test_mols = molecules4
fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints
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
        assert Fingerprint(property="Contains C",
                    noun="molecule",
                    verb="are",
                    label_type="binary",
                    calculator="add",
                    mol_format="smiles"), "Fingerprint generation did not work properly"

"""class TestMolSet(unittest.TestCase):
    global analysis
    def test_molset_generation(self):
            analysis = MolSet(molecules4,
                            fingerprints = fingerprints,
                            mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                            significance=0.1,
                            file="exmoset/data/QM9_Data.csv")"""


class TestGetOutliers(unittest.TestCase):
    def test_outlier_identification(self):
        molspace = MolSpace(molecules=molecules2,
                        fingerprints = fingerprints,
                        mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                        significance=0.1)
        print(molspace.calc_vector(molspace["Full"]))
        print(molspace.get_outliers(molspace.clusters["Full"]))
        assert molspace.get_outliers(molspace.clusters["Full"]) == np.array([5]), "Outlier identification is not working correctly."


class TestMolSpace(unittest.TestCase):
    def test_molspace_gen(self):
        space = MolSpace(fingerprints = fingerprints,
                         molecules=molecules3,
                         mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                         significance=0.1,
                         file="exmoset/data/QM9_Data.csv",
                         index_col="SMILES")

class TestMolFileSpace(unittest.TestCase):
    global filemolspace
    filemolspace = MolSpace(fingerprints=fingerprints,
                         file="exmoset/data/QM9_tiny.csv",
                         mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                         index_col="SMILES",
                         clusters={"Test" : np.array([0,0,0,0,1,1,1,1,1])})
    def test_molspace_clusters(self):
        assert any(filemolspace.clusters["Test"][0]), "Clustering is not working"

if __name__ == "__main__":
    unittest.main()
