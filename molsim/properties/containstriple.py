import numpy as np
from rdkit import Chem
from .property import Property

class ContainsTriple(Property):
    """
    How many rings the molecule contains
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)

        self.values = self.calc_property(molecules)

    @staticmethod
    def maxbondtype(mol):
        types = [b.GetBondType() for b in mol.GetBonds()]
        if "TRIPLE" in types:
            return 3
        elif "DOUBLE" in types:
            return 2
        else:
            return 1

    def calc_property(self,molecules):
        return np.array([self.maxbondtype(x) for x in molecules])

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            return f"Contains Triple Bond"
