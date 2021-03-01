"""
Abstract Class for individual molecular properties
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np
from math import log

class Property(ABC):
    """
    Abstract Class that handles individual molecular properties to provide extensability and
    customize print operations
    """
    def __init__(self):
        pass

    def __str__(self):
        return f"{np.array_repr(self.values)}"

    def name(self):
        return type(self).__name__

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

    def entropy(self,base=None):
        """ Computes entropy of label distribution. """

        n_labels = len(self.values)

        if n_labels <= 1:
            return 0

        value,counts = np.unique(self.values, return_counts=True)
        probs = counts / n_labels
        n_classes = np.count_nonzero(probs)
        if n_classes <= 1:
            return 0
        ent = 0.
        # Compute entropy
        base = 2 if base is None else base
        for i in probs:
            ent -= i * log(i, base)

        return ent

    @abstractmethod
    def summative_label(self):
        pass

    @abstractmethod
    def calc_property(self):
        pass

    @staticmethod
    def convert_mols(molecules):
        return [Chem.MolFromSmiles(mol) for mol in molecules]
