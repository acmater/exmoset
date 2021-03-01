import numpy as np
from rdkit import Chem
from property import Property


class BondType(Property):
    """
    Whether or not the molecule contains Nitrogen
    """
    def __init__(self,bond_type,molecules):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)
        self.bond_type   = bond_type
        self.values      = self.calc_property(molecules,bond_type)

    @staticmethod
    def calc_property(molecules,bond_type):
        return np.array([1 if bond_type in [b.GetBondType().__str__() for b in mol.GetBonds()] else 0 for mol in molecules])

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            return f"Contains {self.bond_type} Bonds" if np.mean(self.values > 0.5) else f"Doesn't Contain {self.bond_type} Bonds"

# Instances

class ContainsSingle(BondType):
    def __init__(self,molecules,bond_type="SINGLE"):
        super().__init__(bond_type,molecules)
        self.values = super().calc_property(molecules,bond_type)

class ContainsDouble(BondType):
    def __init__(self,molecules,bond_type="DOUBLE"):
        super().__init__(bond_type,molecules)
        self.values = super().calc_property(molecules,bond_type)

class ContainsTriple(BondType):
    def __init__(self,molecules,bond_type="TRIPLE"):
        super().__init__(bond_type,molecules)
        self.values = super().calc_property(molecules,bond_type)
