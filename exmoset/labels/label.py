"""
Abstract Class for molecular labels
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np
from math import log, e
from entropy_estimators import continuous

class Label(ABC):
    def __init__(self):
        pass

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

    @staticmethod
    def neg(string):
        if "are" in string:
            return string.replace("are", "are not")
        elif "contain" in string:
            return string.replace("contain", "do not contain")

    def entropy(self,base=None):
        """ Computes entropy of label distribution.

        Parameters
        ----------
        values : np.ndarray
            The array of label values.

        base: float
            Floating number used as base for entropy calculation.
        """
        if type(self).__name__ == "Continuous":
            ent = continuous.get_h(self.values,k=10,norm="euclidean",min_dist=0.001)

        else:
            n_labels = len(self.values)

            if n_labels <= 1:
                return 0

            value, counts = np.unique(self.values, return_counts=True)
            probs = counts / n_labels
            n_classes = np.count_nonzero(probs)
            if n_classes <= 1:
                return 0
            ent = 0.
            # Compute entropy
            base = e if base is None else base
            for i in probs:
                ent -= i * log(i, base)

        return ent

    @abstractmethod
    def summary(self):
        pass
