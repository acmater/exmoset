import numpy as np
from rdkit import Chem
from property import Property

class ZPVE(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules,df):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.ent_type    = "Continuous"
        self.values      = self.calc_property(molecules,df)

    @staticmethod
    def calc_property(molecules,df):
        return df["Dipole Moment"].to_numpy()

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            print(self.entropy())
            return f"Aromatic" if np.mean(self.values > 0.5) else f"Non-Aromatic"
