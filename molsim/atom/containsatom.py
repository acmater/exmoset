import numpy as np
from rdkit import Chem
from molsim.property import Property

class ContainsAtom(Property):
    """
    Whether or not the molecule contains Nitrogen
    """
    def __init__(self,molecules,atoms,df=None):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)
        self.atoms    = atoms
        self.values   = self.calc_property(molecules,atoms)
        self.ent_type = "Discrete"

    def calc_property(self,molecules,atoms):
        atoms = {}
        for atom in self.atoms:
            atoms[atom] = np.array([1 if atom in [a.GetSymbol() for a in mol.GetAtoms()] else 0 for mol in molecules])
        return atoms

    def summative_label(self,significance=0.1,verbose=True):
        summary = []
        for atom in self.atoms:
            if self.entropy(self.values[atom],self.ent_type) < significance:
                if verbose:
                    print(f"Working on Atom: {atom}")
                    print(f"Inside the inner loop. Entropy is {self.entropy(self.values[atom],self.ent_type)}")
                    print(f"Due to signifiance, calculating average presence {atom}")
                    print(f"Average value of {atom} is {np.mean(self.values[atom])}")
                    print()

                    summary.append(f"Contains {atom}" if np.mean(self.values[atom]) > 0.5 else f"Does not contain {atom}")

        if summary:
            return "\n".join(summary)
