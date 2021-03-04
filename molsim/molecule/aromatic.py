import numpy as np
from rdkit import Chem
from molsim.property import Property

class Aromatic(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules,df=None):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values  = self.calc_property(molecules)
        self.ent_type    = "Discrete"

    @staticmethod
    def calc_property(molecules):
        return np.array([bool(x.GetAromaticAtoms()) for x in molecules],dtype=np.int)

    def summative_label(self,significance=0.1,verbose=False):
        if self.entropy(self.values) < significance:
            return f"Aromatic" if np.mean(self.values > 0.5) else f"Non-Aromatic"
