"""
Abstract Class for individual molecular properties
"""
from abc import ABC, abstractmethod
from rdkit import Chem
import numpy as np
from math import log, e
from entropy_estimators import continuous

class Property(ABC):
    """
    Abstract Class that handles individual molecular properties to provide extensability and
    customize print operations
    """
    def __init__(self,molecules):
        assert isinstance(molecules[0],str) or isinstance(molecules[0],Chem.rdchem.Mol), "Not a supported molecular format"

    def __str__(self):
        return f"{np.array_repr(self.values)}"

    def name(self):
        """
        Returns name of property subclass.
        """
        return type(self).__name__

    def __getitem__(self,idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

    def entropy(self,values,base=None):
        """ Computes entropy of label distribution.

        Parameters
        ----------
        values : np.ndarray
            The array of label values.

        base: float
            Floating number used as base for entropy calculation.
        """
        if self.ent_type == "Discrete":
            n_labels = len(values)

            if n_labels <= 1:
                return 0

            value,counts = np.unique(values, return_counts=True)
            probs = counts / n_labels
            n_classes = np.count_nonzero(probs)
            if n_classes <= 1:
                return 0
            ent = 0.
            # Compute entropy
            base = e if base is None else base
            for i in probs:
                ent -= i * log(i, base)

        elif self.ent_type == "Continuous":
            ent = continuous.get_h(values,k=10,norm="euclidean",min_dist=0.001)

        return ent

    @staticmethod
    def convert_mols(molecules,debug=False):
        """
        Helper function to convert smiles into Chem.rdchem.Mol objects

        Parameters
        ----------
        molecules : [str]
            Iterable of SMILES representations of the molecules of interest.

        debug : bool, default=False
            Whether or not failed smiles will also be returned for debugging purposes.
        """

        failed     = []
        successful = []
        for mol in molecules:
            molobj = Chem.MolFromSmiles(mol)
            if molobj is not None:
                successful.append(molobj)
            else:
                failed.append(molobj)
                print(f"Failed to convert {mol}")
        if debug:
            return successful, failed
        else:
            return successful

    def summative_label(self):
        # I think I can make this code general. Just need to determine a few things.
        summary = {}
        for prop in self.atoms:
            if self.entropy(self[atom]) < significance:
                if verbose:
                    print(f"Working on Atom: {atom}")
                    print(f"Inside the inner loop. Entropy is {self.entropy(self[atom])}")
                    print(f"Due to signifiance, calculating average presence {atom}")
                    print(f"Average value of {atom} is {np.mean(self[atom])}")
                    print()

                summary.append(f"Contains {atom}" if np.mean(self[atom]) > 0.5 else f"Does not contain {atom}")

        if summary:
            return "\n".join(summary)


    @abstractmethod
    def calc_property(self):
        """
        Method to calculate the label associated with each property for each molecule within the subset.
        """
        pass
