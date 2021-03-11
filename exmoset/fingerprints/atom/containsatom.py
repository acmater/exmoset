import numpy as np
from ..fingerprint import Fingerprint

atoms = ["C","N","O","F"]

def contains(atom):
    def sub_contains(mol):
        return 1 if atom in mol else 0
    return sub_contains

atom_fingerprints = []

for atom in atoms:
    atom_fingerprints.append(Fingerprint(name=f"Contains {atom}",
                context="Molecules",
                label_type="binary",
                calculator=contains(atom),
                mol_format="smiles"))
