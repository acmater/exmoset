import numpy as np
from rdkit import Chem
from property import Property

class Aromatic(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"
        if isinstance(molecules[0],str):
            print("Using")
            self.molecules = super().convert_mols(molecules)
        else:
            self.molecules = molecules

        self.values = self.calc_property()

    def calc_property(self):
        return [bool(x.GetAromaticAtoms()) for x in self.molecules]

    def __str__(self):
        return f"{self.molecules}"

    def __getitem__(self,idx):
        return self.molecules[idx]



if __name__ == "__main__":
    a = "CCNC"
    b = "c1ccccc1"
    mols = list(map(Chem.MolFromSmiles,[a,b]))
    print(Aromatic(mols).values)
