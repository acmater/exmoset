import numpy as np
from rdkit import Chem
from molsim.property import Property

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

    def summative_label(self,significance=0.1,verbose=False):
        summary = []
        for substructure in self.substructures:
            if verbose:
                print(f"Working on structure: {substructure}")
                if self.entropy(selfs[substructure]) < significance:
                    print(f"Inside the inner loop. Entropy is {self.entropy(self[substructure])}")
                    print(f"Due to signifiance, identifying presence of absence of {substructure}")
                    print(f"Average value of {substructure} is {np.mean(self[substructure])}")
                    print()
                summary.append(f"Contains {substructure} Group" if np.mean(self[substructure]) > 0.5 else f"Doesn't Contain {substructure} Group")


        if summary:
            return "\n".join(summary)
