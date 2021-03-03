import numpy as np
from rdkit import Chem
import pandas as pd

from molecule import *
from atom import ContainsAtom
from bond import ContainsBond
from substructure import Substructure
from data import *

properties = [Aromatic,NumRings,NumAtoms]

class Similarity_Analysis():
    def __init__(self,molecules,
                      properties,
                      atoms=None,
                      bonds=None,
                      molprops=None,
                      substructures=None,
                      significance=0.1,
                      file=None):
        self.molecules = [Chem.MolFromSmiles(mol) for mol in molecules]
        if file is not None:
            print(f"Importing {file}")
            df     = pd.read_csv(file,index_col="SMILES")
            sub_df = df.loc[molecules]

        if atoms is not None:
            atom_props = ContainsAtom(molecules,atoms)

        if bonds is not None:
            bond_props = ContainsBond(molecules,bonds)

        if molprops is not None:
            molecular_properties = MolProp(molecules,molprops,sub_df)

        if substructures is not None:
            substructures = Substructure(molecules,substructures)



        self.properties = [prop(self.molecules,df=sub_df) for prop in properties] + [atom_props,bond_props,molecular_properties,substructures]
        self.significance = significance

    def __str__(self):
        header  = ["Subset Description"]
        labels = [prop.summative_label(significance=self.significance) for prop in self.properties]
        results = ["\t{:<60}".format(label) for label in labels if label is not None]

        if len(results) == 0:
            return "No Discernible pattern"
        else:
            final = "\n".join(header + results)
            return final


if __name__ == "__main__":
    analy = Similarity_Analysis(molecules6,
                                properties,
                                atoms=["C","N","O","F"],
                                bonds=["SINGLE","DOUBLE","TRIPLE"],
                                molprops=["Dipole Moment","Isotropic Polarizability", "Electronic Spatial Extent"],
                                substructures = ["[OH]","[NH2]","[CC]"],
                                significance=0.1,
                                file="data/QM9_Data.csv")
    print(analy)

    df = pd.read_csv("data/QM9_Data.csv",index_col="SMILES")
    x = df.loc[molecules6]["Dipole Moment"].to_numpy()
    import matplotlib.pyplot as plt
    plt.hist(x)
    plt.show()
