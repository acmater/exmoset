import numpy as np
from rdkit import Chem
from .property import Property

class Aromatic(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values = self.calc_property(molecules)

    @staticmethod
    def calc_property(molecules):
        return np.array([bool(x.GetAromaticAtoms()) for x in molecules],dtype=np.int)

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            return f"Aromatic" if np.mean(self.values > 0.5) else f"Non-Aromatic"


if __name__ == "__main__":
    a = "CCNC"
    b = "c1ccccc1"
    mols = list(map(Chem.MolFromSmiles,[a,b]))
    print(Aromatic(mols).values)
