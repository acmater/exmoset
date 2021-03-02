import numpy as np
from rdkit import Chem
from property import Property

class Substructure(Property):
    def __init__(self,molecules,substructures):
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
        self.substructures  = substructures
        self.values         = self.calc_property(molecules,substructures)
        self.ent_type       = "Discrete"

    @staticmethod
    def calc_property(molecules,substructures):
        matches = {}
        for substructure in substructures:
            substructure_patt = Chem.MolFromSmarts(substructure)
            matches[substructure] = np.array([mol.HasSubstructMatch(substructure_patt) for mol in molecules],dtype=np.int)
        return matches

    def summative_label(self,significance=0.1):
        for substructure in self.substructures:
            if self.entropy() < significance:
                # What in the world is going on here?!?!?!
                print(np.mean(self.values[substructure]) > 0.5)
                if np.mean(self.values[substructure]) > 0.5 == False:
                    print("This Value is False")
                print(np.mean(self.values[substructure]) > 0.5)
                return f"Contains {substructure} Group" if np.mean(self.values[substructure] > 0.5) else f"Doesn't Contain {substructure} Group"
