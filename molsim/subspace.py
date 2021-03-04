import numpy as np
from rdkit import Chem
import pandas as pd

from molecule import *
from atom import ContainsAtom
from bond import ContainsBond
from substructure import Substructure
from data import *

properties = [Aromatic]#,NumRings,NumAtoms]

class Subspace():
    def __init__(self,molecules,
                      properties,
                      atoms=None,
                      bonds=None,
                      molprops=None,
                      substructures=None,
                      verbose=False,
                      significance=0.1,
                      file=None):
        self.molecules, self.smiles = self.convert_mols(molecules)
        self.verbose = verbose

        if file is not None:
            print(f"Importing {file}")
            df     = pd.read_csv(file,index_col="SMILES")
            sub_df = df.loc[self.smiles]

        if atoms is not None:
            atom_props = ContainsAtom(self.molecules,atoms)

        if bonds is not None:
            bond_props = ContainsBond(self.molecules,bonds)

        if molprops is not None:
            molecular_properties = MolProp(self.smiles,molprops,sub_df)

        if substructures is not None:
            substructures = Substructure(self.molecules,substructures)

        self.properties = [prop(self.molecules,df=sub_df) for prop in properties] + [atom_props,bond_props,molecular_properties,substructures]
        self.significance = significance



    def __str__(self):
        header  = ["Subset Description"]
        labels = [prop.summative_label(significance=self.significance,verbose=self.verbose) for prop in self.properties]
        results = ["\t{:<60}".format(label) for label in labels if label is not None]

        if len(results) == 0:
            return "No Discernible pattern"
        else:
            final = "\n".join(header + results)
            return final

    def __or__(self, other):
        pass

    @staticmethod
    def convert_mols(molecules,debug=False):
        failed     = []
        successful = []
        smiles     = []
        for mol in molecules:
            molobj = Chem.MolFromSmiles(mol)
            if molobj is not None:
                successful.append(molobj)
                smiles.append(mol)
            else:
                failed.append(molobj)
                print(f"Failed to convert {mol}")
        if debug:
            return successful, smiles, failed
        else:
            return successful, smiles


if __name__ == "__main__":
    analy = Subspace(molecules10,
                        properties,
                        atoms=["C","N","O","F"],
                        bonds=["SINGLE","DOUBLE","TRIPLE"],
                        molprops=["Dipole Moment","Isotropic Polarizability", "Electronic Spatial Extent", "Rotational Constant A"],
                        substructures = ["[OH]","[NH2]","[CC]"],
                        significance=0.1,
                        file="data/QM9_Data.csv")
    print(analy)
