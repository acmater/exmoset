import numpy as np
from exmoset.property import Property
from exmoset.labels import Binary

class ContainsAtom(Property):
    """
    Whether or not the molecule contains Nitrogen
    """
    def __init__(self,molecules,atoms):
        """
        Initialization of class

        Paramters
        ---------
        molecules : [str] or [Chem.rdchem.Mol]
            A list of molecules to be analyzed.

        atoms : [str]
            A list of strings representing atom types to be considered.

        """
        super().__init__(molecules)
        if isinstance(molecules[0],str):
            molecules = super().convert_mols(molecules)
        self.atoms    = atoms
        self.values   = self.calc_property(molecules)
        self.ent_type = "Discrete"

    def calc_property(self,molecules):
        """
        Calculates the presence or absence of each atom type in the set of molecules provided.

        Returns
        -------
        atoms : {"<atom_type>" : np.array}
            Each atom type (expressed as elemental symbol strings is mapped to the binary vector describing
            whether or not it is present in each molecule.
        """
        atoms = {}
        for atom in self.atoms:
            atoms[atom] = Binary(atom,
                                 np.array([1 if atom in [a.GetSymbol() for a in mol.GetAtoms()] else 0 for mol in molecules]),
                                 context="part")
        return atoms

    def summative_label(self,significance=0.1,verbose=True):
        summary = []
        for atom, label in self.values.items():
            if self.entropy(self[atom]) < significance:
                if verbose:
                    print(f"Working on Atom: {atom}")
                    print(f"Inside the inner loop. Entropy is {self.entropy(self[atom])}")
                    print(f"Due to signifiance, calculating average presence {atom}")
                    print(f"Average value of {atom} is {np.mean(self[atom])}")
                    print()

                summary.append(label.summary())

        if summary:
            return "\n".join(summary)
