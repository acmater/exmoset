"""
Abstract Class for individual molecular properties
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np

class Property(ABC):
    """
    Abstract Class that handles individual molecular properties to provide extensability and
    customize print operations
    """
    def __init__(self):
        pass

    def __str__(self):
        return f"{np.array_repr(self.values)}"

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

    @abstractmethod
    def calc_property(self):
        pass

    @staticmethod
    def convert_mols(molecules):
        return [Chem.MolFromSmiles(mol) for mol in molecules]
