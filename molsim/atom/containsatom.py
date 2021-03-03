import numpy as np
from rdkit import Chem
from molsim.property import Property

print("worked this time")
print()
class ContainsAtom(Property):
    """
    Whether or not the molecule contains Nitrogen
    """
    def __init__(self,atom,molecules,df=None):
        super().__init__(molecules)
        self.atom   = atom
        self.values = self.calc_property(molecules,atom)
        self.ent_type = "Discrete"

    @staticmethod
    def calc_property(molecules,atom):
        return np.array([1 if atom in [a.GetSymbol() for a in mol.GetAtoms()] else 0 for mol in molecules])

    def summative_label(self,significance=0.1):
        if self.entropy(self.values,self.ent_type) < significance:
            return f"Contains {self.atom}" if np.mean(self.values) > 0.5 else f"Doesn't Contain {self.atom}"

# Instances

class ContainsNitrogen(ContainsAtom):
    def __init__(self,molecules,atom="N",df=None):
        super().__init__(atom,molecules)
        self.values = super().calc_property(molecules,atom)

class ContainsCarbon(ContainsAtom):
    def __init__(self,molecules,atom="C",df=None):
        super().__init__(atom,molecules)
        self.values = super().calc_property(molecules,atom)

class ContainsOxygen(ContainsAtom):
    def __init__(self,molecules,atom="O",df=None):
        super().__init__(atom,molecules)
        self.values = super().calc_property(molecules,atom)

class ContainsFluorine(ContainsAtom):
    def __init__(self,molecules,atom="F",df=None):
        super().__init__(atom,molecules)
        self.values = super().calc_property(molecules,atom)
