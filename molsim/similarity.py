import numpy as np
from rdkit import Chem
import pandas as pd

from molecule import *
from atom import *
from bond import *
from substructure import *
from data import *

properties = [Aromatic,NumRings,NumAtoms,
ContainsNitrogen,
ContainsCarbon,
ContainsFluorine,
ContainsOxygen,
ContainsSingle,
ContainsDouble,
ContainsTriple]

class Similarity_Analysis():
    def __init__(self,molecules,
                      properties,
                      molprops=None,
                      substructures=None,
                      significance=0.1,
                      file=None):
        self.molecules = [Chem.MolFromSmiles(mol) for mol in molecules]
        if file is not None:
            print(f"Importing {file}")
            df     = pd.read_csv(file,index_col="SMILES")
            sub_df = df.loc[molecules]

        if molprops is not None:
            molecular_properties = MolProp(molecules,molprops,sub_df)

        if substructures is not None:
            substructures = Substructure(molecules,substructures)

        self.properties = [prop(self.molecules,df=sub_df) for prop in properties] + [substructures,molecular_properties]
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
    analy = Similarity_Analysis(molecules10,
                                properties,
                                molprops=["Dipole Moment","Isotropic Polarizability"],
                                substructures = ["[OH]","[NH2]","[CC]"],
                                significance=0.1,
                                file="data/QM9_Data.csv")
    print(analy)
