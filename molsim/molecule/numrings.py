import numpy as np
from rdkit import Chem
from property import Property

class NumRings(Property):
    """
    How many rings the molecule contains
    """
    def __init__(self,molecules,df=None):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values = self.calc_property(molecules)
        self.ent_type = "Discrete"

    @staticmethod
    def calc_property(molecules):
        return np.array([x.GetRingInfo().NumRings() for x in molecules])

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            return f"{int(np.round(np.mean(self.values)))} Rings"
