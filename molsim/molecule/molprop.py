import numpy as np
from rdkit import Chem
from property import Property

class MolProp(Property):
    """
    Whether or not a molecular is aromatic
    """
    def __init__(self,molecules,properties,df):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.ent_type    = "Continuous"
        self.properties  = properties
        self.values      = self.calc_property(molecules,df)

    @staticmethod
    def calc_property(molecules,properties,df):
        props = {}
        for prop in properties:
            props[prop] = df[prop].to_numpy()
        return props

    def summative_label(self,significance=0.1):
        for prop in properties:
            if self.entropy() < significance:
                print(self.entropy())
                return f"Aromatic" if np.mean(self.values > 0.5) else f"Non-Aromatic"
