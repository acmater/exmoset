import numpy as np
from molsim.property import Property

class ContainsBond(Property):
    """
    Whether or not the molecule contains Nitrogen
    """
    def __init__(self,molecules,bonds):
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)
        self.bonds    = bonds
        self.values   = self.calc_property(molecules)
        self.ent_type = "Discrete"

    def calc_property(self,molecules):
        bonds = {}
        # This really slow, change how the code functions to generate the bond description and then test each type against it.
        for bond in self.bonds:
            bonds[bond] = np.array([1 if bond in [b.GetBondType().__str__() for b in mol.GetBonds()] else 0 for mol in molecules])
        return bonds

    def summative_label(self,significance=0.1,verbose=True):
        summary = []
        for bond in self.bonds:
            if self.entropy(self[bond]) < significance:
                if verbose:
                    print(f"Working on Bond: {bond}")
                    print(f"Inside the inner loop. Entropy is {self.entropy(self[bond])}")
                    print(f"Due to signifiance, calculating average presence {bond}")
                    print(f"Average value of {bond} is {np.mean(self[bond])}")
                    print()

                    summary.append(f"Contains {bond}" if np.mean(self[bond]) > 0.5 else f"Does not contain {bond}")

        if summary:
            return "\n".join(summary)
