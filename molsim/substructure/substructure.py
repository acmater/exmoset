import numpy as np
from rdkit import Chem
from property import Property

class Substructure(Property):
    def __init__(self,molecules,substructure):
        """
        Initializes the parameters of the class

        Parameters
        ----------
        molecules : np.array[str,Chem.rdchem.Mol]

            The iterable of molecules that should be examined

        substructre : rdkit.Chem.rdchem.Mol

            A valid SMARTS string that will be queried against each molecule in molecules
        """
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)
        self.substructure  = substructure
        self.values        = self.calc_property(molecules,substructure)

    @staticmethod
    def calc_property(molecules,substructure):
        substructure_patt = Chem.MolFromSmarts(substructure)
        return np.array([mol.HasSubstructMatch(substructure_patt) for mol in molecules],dtype=np.int)

    def summative_label(self,significance=0.1):
        if self.entropy() < significance:
            # What in the world is going on here?!?!?!
            print(np.mean(self.values) > 0.5)
            if np.mean(self.values) > 0.5 == False:
                print("This Value is False")
            print(np.mean(self.values) > 0.5)
            return f"Contains {self.substructure} Group" if np.mean(self.values > 0.5) else f"Doesn't Contain {self.substructure} Group"

# Instances

class ContainsAlcohol(Substructure):
    def __init__(self,molecules,substructure="[OH]"):
        super().__init__(molecules,substructure)
        self.values = super().calc_property(molecules,substructure)

class ContainsPrimaryAmine(Substructure):
    def __init__(self,molecules,substructure="[NH2]"):
        super().__init__(molecules,substructure)
        self.values = super().calc_property(molecules,substructure)
