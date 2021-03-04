"""
Abstract Class for molecular labels
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np

class Label(ABC):
    def __init__(self):
        pass

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)
