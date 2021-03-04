import numpy as np
from rdkit import Chem
from exmoset.abstract import Property
from exmoset.labels import Binary

class Aromatic(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules,df=None):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.label   = self.calc_property(molecules)

    @staticmethod
    def calc_property(molecules):
        return Binary("Aromatic",
                      np.array([bool(x.GetAromaticAtoms()) for x in molecules],dtype=np.int),
                      context="whole")

    def summative_label(self,significance=0.1,verbose=False):
        if self.entropy(self.label) < significance:
            return self.label.summary()
