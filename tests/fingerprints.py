"""
Prewritten fingerprints for atom based properties.
"""

import numpy as np
from exmoset.fingerprints import Fingerprint
from rdkit import Chem

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

"""
Prewritten fingerprints for bond based properties.
"""

bonds = ["SINGLE","DOUBLE","TRIPLE"]

def contains(bond):
    def sub_contains(mol):
         return 1 if bond in [b.GetBondType().__str__() for b in mol.GetBonds()] else 0
    return sub_contains

bond_fingerprints = []

for bond in bonds:
    bond_fingerprints.append(Fingerprint(name=f"Contains {bond}",
                context="Molecules",
                label_type="binary",
                calculator=contains(bond),
                mol_format="rd"))

"""
General molecular fingerprints
"""

def aromatic(mol):
    return np.array([bool(mol.GetAromaticAtoms())],dtype=np.int)
def num_atoms(mol):
    return np.array([len(mol.GetAtoms())])
def num_rings(mol):
    return np.array([mol.GetRingInfo().NumRings()])


general_fingerprints = [Fingerprint(name="Aromatic",
                                context="whole",
                                label_type="binary",
                                calculator=aromatic,
                                mol_format="rd"),

                    Fingerprint(name="Number of Atoms",
                                context="whole",
                                label_type="multiclass",
                                calculator=num_atoms,
                                mol_format="rd"),

                    Fingerprint(name="Number of Rings",
                                context="whole",
                                label_type="multiclass",
                                calculator=num_rings,
                                mol_format="rd")]


"""
Prewritten fingerprints for substructure matching based properties.
"""

substructures = ["[CC]","[OH]","[NH2]"]

def contains(substructure):
    def sub_contains(mol):
        return np.array([mol.HasSubstructMatch(Chem.MolFromSmarts(substructure))],dtype=np.int)
    return sub_contains

substructure_fingerprints = []

for substructure in substructures:
    substructure_fingerprints.append(Fingerprint(name=f"Contains {substructure}",
                context="Molecules",
                label_type="binary",
                calculator=contains(substructure),
                mol_format="rd"))
