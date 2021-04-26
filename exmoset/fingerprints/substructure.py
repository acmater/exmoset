"""
Prewritten fingerprints for substructure matching based properties.
"""

import numpy as np
from .fingerprint import Fingerprint
from rdkit import Chem

substructures = ["[CC]","[OH]","[NH2]","*=O","*O*"]

def contains(substructure):
    match = Chem.MolFromSmarts(substructure)
    def sub_contains(mol):
        return int(mol.HasSubstructMatch(Chem.MolFromSmarts(match)))
    return sub_contains

substructure_fingerprints = []

for substructure in substructures:
    substructure_fingerprints.append(Fingerprint(property=f"{substructure}",
                verb="contain",
                noun="Molecules",
                label_type="binary",
                calculator=contains(substructure),
                mol_format="rd"))
