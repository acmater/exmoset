"""
Abstract Class for individual molecular properties
"""
from abc import ABC, abstractmethod
from rdkit import Chem

class Property(ABC):
    """
    Abstract Class that handles individual molecular properties to provide extensability and
    customize print operations
    """
    def __init__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def calc_property(self):
        pass

    @staticmethod
    def convert_mols(molecules):
        return [Chem.MolFromSmiles(mol) for mol in molecules]
