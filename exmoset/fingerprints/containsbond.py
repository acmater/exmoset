"""
Prewritten fingerprints for bond based properties.
"""

import numpy as np
from .fingerprint import Fingerprint

bonds = ["SINGLE","DOUBLE","TRIPLE"]

def contains(bond):
    def sub_contains(mol):
         return 1 if bond in [b.GetBondType().__str__() for b in mol.GetBonds()] else 0
    return sub_contains

bond_fingerprints = []

for bond in bonds:
    bond_fingerprints.append(Fingerprint(property=f"{bond.lower()} bonds",
                verb="contain",
                context="Molecules",
                label_type="binary",
                calculator=contains(bond),
                mol_format="rd"))
