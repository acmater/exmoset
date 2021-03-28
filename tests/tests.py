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

"""
class TestBinaryProperty(unittest.TestCase):
    global binary
    binary = Binary("binary",np.ones((100,1)))
    def test_av(self):
        assert binary.av == 1, "Binary label averaging not working."
    def test_entropy(self):
        assert binary.entropy == 0, "Entropy testing for Binary labels not working."
    def test_property(self):
        assert binary.property == "binary", "Binary property (name) not working properly."
    def test_plot(self):
        fig = binary.plot()
        assert isinstance(fig, matplotlib.figure.Figure), "Binary Figure plotting not working correctly."

class TestMulticlassLabel(unittest.TestCase):
    global multi
    multi = Multiclass("multi",np.arange(10))
    def test_av(self):
        assert multi.av == 4, "Multiclass label averaging not working."
    def test_entropy(self):
        assert multi.entropy != 0, "Entropy testing for Multiclass labels not working."
    def test_property(self):
        assert multi.property == "multi", "Multiclass label property (name) not working properly."
    def test_plot(self):
        fig = multi.plot()
        assert isinstance(fig, matplotlib.figure.Figure), "Multiclass Figure plotting not working correctly."

class TestContinuousclassLabel(unittest.TestCase):
    global cont
    cont = Continuous("cont",np.arange(10))
    def test_av(self):
        assert cont.av == 4.5, "Continuous label averaging not working."
    def test_entropy(self):
        assert cont.entropy != 0, "Entropy testing for Continuous labels not working."
    def test_property(self):
        assert cont.property == "cont", "Continuous label property not working properly."
    def test_plot(self):
        fig = cont.plot()
        assert isinstance(fig, matplotlib.figure.Figure), "Continuous Figure plotting not working correctly."

"""
"""class TestMolSet(unittest.TestCase):
    global analysis
    def test_molset_generation(self):
            analysis = MolSet(molecules4,
                            fingerprints = fingerprints,
                            mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                            significance=0.1,
                            file="exmoset/data/QM9_Data.csv")"""

class TestMolFileSpace(unittest.TestCase):
    global filemolspace
    filemolspace = MolSpace(fingerprints=fingerprints,
                         file="exmoset/data/QM9_tiny.csv",
                         mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                         index_col="SMILES",
                         clusters={"Test" : np.array([0,0,0,0,1,1,1,1,1])})
    def test_molspace_clusters(self):
        assert any(filemolspace.clusters["Test"][0]), "Clustering is not working"


class TestGetOutliers(unittest.TestCase):
    def test_outlier_identification(self):
        molspace = MolSpace(molecules=molecules2,
                        fingerprints = fingerprints,
                        mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                        significance=0.1)
        assert molspace.get_outliers(molspace.clusters["Full"]) == np.array([5]), "Outlier identification is not working correctly."


class TestMolSpace(unittest.TestCase):
    def test_molspace_gen(self):
        space = MolSpace(fingerprints = fingerprints,
                         molecules=molecules3,
                         mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                         significance=0.1,
                         index_col="SMILES")

class TestSpaceUpdating(unittest.TestCase):
    def test_update_fp(self):
        def conjugated(mol):
            return 1 if any([b.GetIsConjugated() for b in mol.GetBonds()]) else 0
        fp = Fingerprint(property="Conjugated",noun="Molecules",verb="are",label_type="binary",calculator=conjugated,mol_format="rd")
        filemolspace.add_fingerprint(fp)
        assert "Conjugated" in filemolspace.fingerprints.keys(), "Fingerprint no longer updating correctly."
    def test_update_cluster(self):
        """
        Tests the update cluster function.
        """
        filemolspace.add_cluster(("Random",np.random.randint(5,size=(3,))))
        assert "Random" in filemolspace.clusters, "Cluster updating not working correctly."

class TestPlottingFunctionality(unittest.TestCase):
    def test_binary_property_plot(self):
        fig = filemolspace.plot("Test", "Aromatic")
        ax = fig.gca()
        p = ax.patches
        assert p[1].get_height() == 4, "Data not plotting correctly for a binary label."
        plt.close()

    def test_multiclass_property_plot(self):
        fig = filemolspace.plot("Test", "Atoms")
        ax = fig.gca()
        p = ax.patches
        assert p[0].get_height() == 3, "Data not plotting correctly for a multiclass label."
        plt.close()

        # TODO Finish
    #def test_continuous_property_plot(self):
    #    fig = filemolspace.plot(filemolspace["Test"][0], "Dipole Moment")
        #plt.show()
    #    print(dir(fig))
        #p = ax.patches
        #assert p[0].get_height() == 3, "Data not plotting correctly for a multiclass label."
        #plt.close()




if __name__ == "__main__":
    unittest.main()
