"""
Prewritten fingerprints for substructure matching based properties.
"""

import numpy as np
from .fingerprint import Fingerprint
from rdkit import Chem

substructures = ["[CC]","[OH]","[NH2]"]

def contains(substructure):
    def sub_contains(mol):
        return np.array([mol.HasSubstructMatch(Chem.MolFromSmarts(substructure))],dtype=np.int)
    return sub_contains

substructure_fingerprints = []

for substructure in substructures:
    substructure_fingerprints.append(Fingerprint(property=f"{substructure}",
                verb="contain",
                noun="Molecules",
                label_type="binary",
                calculator=contains(substructure),
                mol_format="rd"))
