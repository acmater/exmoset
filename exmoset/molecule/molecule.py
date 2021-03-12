import numpy as np
from rdkit import Chem

class Molecule():
    """
    Python abstract class for a molecule. Protein is
    initialized with a variety of properties and utilises dictionary type syntaxing
    to make it compliant with the remainder of the code.

    Parameters
    ----------
    rep :

        The only required argument. Is used to calculate length.
    """
    def __init__(self,rep,**kwargs):
        self.rep = rep
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        parsed = []
        for key, value in vars(self).items():
            if isinstance(value,np.ndarray):
                parsed.append(f"{key}=np.array({list(value)})")
            elif isinstance(value, str):
                parsed.append(f"{key}='{value}'")
            else:
                parsed.append(f"{key}={value}")
        return f"Molecule(" + ",".join(parsed) + ")"

    def __getitem__(self,keys):
        return self.__dict__[keys]

    def __len__(self):
        return len(self.rep)

    def __eq__(self,other):
        return True if self.rep == other.rep else False

if __name__ == "__main__":
    mol_converters = {"mol" : Chem.MolFromSmiles("CCC")}
    test = Molecule("CCC",**mol_converters)
    test2 = Molecule("CCC")
