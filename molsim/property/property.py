"""
Abstract Class for individual molecular properties
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np
from math import log, e
from entropy_estimators import continuous

class Property(ABC):
    """
    Abstract Class that handles individual molecular properties to provide extensability and
    customize print operations
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"

    def __str__(self):
        return f"{np.array_repr(self.values)}"

    def name(self):
        return type(self).__name__

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

    def entropy(self,values,base=None):
        """ Computes entropy of label distribution. """
        if self.ent_type == "Discrete":
            n_labels = len(values)

            if n_labels <= 1:
                return 0

            value,counts = np.unique(values, return_counts=True)
            probs = counts / n_labels
            n_classes = np.count_nonzero(probs)
            if n_classes <= 1:
                return 0
            ent = 0.
            # Compute entropy
            base = e if base is None else base
            for i in probs:
                ent -= i * log(i, base)

        elif self.ent_type == "Continuous":
            ent = continuous.get_h(values,k=10,norm="euclidean",min_dist=0.001)

        return ent

    @staticmethod
    def convert_mols(molecules):
        return [Chem.MolFromSmiles(mol) for mol in molecules]

    @abstractmethod
    def summative_label(self):
        pass

    @abstractmethod
    def calc_property(self):
        pass
