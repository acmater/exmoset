import numpy as np
from rdkit import Chem
from property import Property

class Aromatic(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values = self.calc_property(molecules)

    @staticmethod
    def calc_property(molecules):
        return np.array([bool(x.GetAromaticAtoms()) for x in molecules],dtype=np.int)

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            return f"Aromatic" if np.mean(self.values > 0.5) else f"Non-Aromatic"
