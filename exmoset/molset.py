import numpy as np
import numpy.ma as ma
from rdkit import Chem
import pandas as pd
import tqdm

from data import *
from utils import Molecule
from labels import Binary, Multiclass, Continuous
from fingerprints import *

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

    significance : float, default = 0.1
        The default signifiance threshold used when calculating whether or not a particular label is significant.

    file : str, default=None
        An optional file (dataframe) that will be imported by pandas and can be accessed by the fingerprints.

    label_tpyes : {str : <Label Class>}, default = {"binary" : Binary, "multiclass" : Multiclass, "continuous" : Continuous}
        A dictionary of possible label types for indexing.
    """
    def __init__(self,molecules,
                      fingerprints,
                      mol_converters={},
                      significance=0.1,
                      file=None,
                      label_types = {"binary"     : Binary,
                                     "multiclass" : Multiclass,
                                     "continuous" : Continuous}):

        self.Molecules = []
        for mol in tqdm.tqdm(molecules):
            formats = {key : mol_converters[key](mol) for key in mol_converters.keys()}
            self.Molecules.append(Molecule(mol, **formats))

        if file is not None:
            print(f"Importing {file}")
            df     = pd.read_csv(file,index_col="SMILES")
            sub_df = df.loc[[mol.smiles for mol in self.Molecules]]

        self.prop_values = np.zeros((len(fingerprints),len(self.Molecules)))

        for i,molecule in enumerate(self.Molecules):
            for j,fp in enumerate(fingerprints):
                if fp.file is not None:
                    self.prop_values[j,i] = fp.calculator(molecule[fp.mol_format],file=sub_df)
                else:
                    self.prop_values[j,i] = fp.calculator(molecule[fp.mol_format])

        self.label_dict = {}
        for i,fp in enumerate(fingerprints):
            self.label_dict[fp.name] = label_types[fp.label_type](fp.name,self.prop_values[i],fp.context)

        self.significance     = significance
        self.properties       = np.array(list(self.label_dict.keys()))
        self.vector,self.mask = self.calc_vector()

    def calc_vector(self):
        """
        Represents the meaningful properties as a vector. Uses numpy's masking feature to only include meaningful
        values within it.
        """
        vector = np.zeros((len(self.label_dict),))
        mask   = np.zeros((len(self.label_dict),))
        for i,prop in enumerate(self.label_dict.values()):
            vector[i] = prop.av
            if prop.entropy < prop.sensitivity:
                mask[i] = 0
            else:
                mask[i] = 1
        return ma.array(vector, mask=mask), mask


    def get_outliers(self):
        # This is clunky af
        full_mask = np.zeros((len(self.mask),len(self)))
        for i in range(len(self)):
            full_mask[:,i] = self.mask
        masked_vals = ma.array(self.prop_values,mask=full_mask)
        return np.where(np.sum((self.vector.reshape(-1,1) - masked_vals),axis=0) != 0) # Need to update distance formulation

    def __len__(self):
        return len(self.Molecules)

    def __str__(self):
        header  = ["Subset Description"]
        labels = [label.summary() for label in self.label_dict.values()]
        results = ["\t{:<60}".format(label) for label in labels if label is not None]

        if len(results) == 0:
            return "No Discernible pattern"
        else:
            final = "\n".join(header + results)
            return final

    def __and__(self, other):
        assert isinstance(other, MolSet)
        return self.properties[ma.where(self.vector == other.vector)]

if __name__ == "__main__":
    fingerprints =  general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints

    analysis = MolSet(molecules2,
                    fingerprints = fingerprints,
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    significance=0.1,
                    file="data/QM9_Data.csv")

    analysis2 = MolSet(molecules4,
                fingerprints = fingerprints,
                mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                significance=0.1,
                file="data/QM9_Data.csv")

    print(analysis.Molecules[5])
