import numpy as np
from rdkit import Chem
import pandas as pd
import tqdm

from data import *
from abstract import Molecule
from labels import Binary, Multiclass, Continuous
from fingerprints import Fingerprint
from fingerprints.atom import atom_fingerprints
from fingerprints.bond import bond_fingerprints

#np.array([1 if "C" in [a.GetSymbol() for a in mol.GetAtoms()]


def contains_OH(mol):
    return np.array([mol.HasSubstructMatch(Chem.MolFromSmarts("[N]"))],dtype=np.int)

def Dipole_Moment(mol,file):
    return file["Dipole Moment"][mol]
def aromatic(mol):
    return np.array([bool(mol.GetAromaticAtoms())],dtype=np.int)
def num_atoms(mol):
    return np.array([len(mol.GetAtoms())])
def num_rings(mol):
    return np.array([mol.GetRingInfo().NumRings()])

fingerprints =  [   Fingerprint(name="Aromatic",
                                context="whole",
                                label_type="binary",
                                calculator=aromatic,
                                mol_format="rd"),

                    Fingerprint(name="Number of Atoms",
                                context="whole",
                                label_type="multiclass",
                                calculator=num_atoms,
                                mol_format="rd"),

                    Fingerprint(name="Number of Rings",
                                context="whole",
                                label_type="multiclass",
                                calculator=num_rings,
                                mol_format="rd"),
# Substructure fingerprints

                     Fingerprint(name="Contains OH",
                                 context="Molecules",
                                 label_type="binary",
                                 calculator=contains_OH,
                                 mol_format="rd"),

# Continuous Dataset Indexing
                    Fingerprint(name="Dipole Moment",
                                context="Molecules",
                                label_type="continuous",
                                calculator=Dipole_Moment,
                                mol_format="smiles",
                                file=True)] + atom_fingerprints + bond_fingerprints

label_types = {"binary"     : Binary,
               "multiclass" : Multiclass,
               "continuous" : Continuous}

class MolSet():
    """
    A molecular set to be analysed.

    Attributes
    ----------
    molecules : list
        A list of molecules that will be parsed into the Molecule class

    properties : list
        A list of property fingerprints that will be calculated for each system

    mol_converters : dict, default = {}
        A list of molecule converters that will be passed as kwargs to the Molecule object

        Example: {"rdkit" : Chem.MolFromSmiles}
        This will then be provided as keyword arguments to Molecule and given the particular mol as an argument.
        Molecule(mol, rdkit=Chem.MolFromSmiles(mol))

    """
    def __init__(self,molecules,
                      fingerprints,
                      mol_converters={},
                      verbose=False,
                      significance=0.1,
                      file=None):

        self.Molecules = []
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            self.Molecules.append(Molecule(mol, **formats))

        if file is not None:
            print(f"Importing {file}")
            df     = pd.read_csv(file,index_col="SMILES")
            sub_df = df.loc[[mol.smiles for mol in self.Molecules]]

        prop_values = np.zeros((len(fingerprints),len(self.Molecules)))

        for i,molecule in enumerate(self.Molecules):
            for j,fp in enumerate(fingerprints):
                if fp.file is not None:
                    prop_values[j,i] = fp.calculator(molecule[fp.mol_format],file=sub_df)
                else:
                    prop_values[j,i] = fp.calculator(molecule[fp.mol_format])

        self.label_dict = {}
        for i,fp in enumerate(fingerprints):
            self.label_dict[fp.name] = label_types[fp.label_type](fp.name,prop_values[i],fp.context)

        print([label.summary() for label in self.label_dict.values()])

        self.significance = significance

    def __str__(self):
        header  = ["Subset Description"]
        labels = [label.summary() for label in self.label_dict.values()]
        results = ["\t{:<60}".format(label) for label in labels if label is not None]

        if len(results) == 0:
            return "No Discernible pattern"
        else:
            final = "\n".join(header + results)
            return final

    def __or__(self, other):
        pass

    def get_outliers(self):
        pass # TODO Implement

if __name__ == "__main__":
    analysis = MolSet(molecules4,
                    fingerprints = fingerprints,
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    significance=0.1,
                    file="data/QM9_Data.csv")
    print(analysis)
