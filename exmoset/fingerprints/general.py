"""
General molecular fingerprints
"""

import numpy as np
from .fingerprint import Fingerprint
from rdkit import Chem

def Dipole_Moment(mol,file):
    return file["Dipole Moment"][mol]
def Electronic_Spatial_Extent(mol, file):
    return file["Electronic Spatial Extent"][mol]

def aromatic(mol):
    return 1 if bool(mol.GetAromaticAtoms()) else 0
def conjugated(mol):
    return 1 if any([b.GetIsConjugated() for b in mol.GetBonds()]) else 0
def num_atoms(mol):
    return len(mol.GetAtoms())
def num_rings(mol):
    return mol.GetRingInfo().NumRings()
def num_branches(mol):
    return mol.count("(")

general_fingerprints = [Fingerprint(property="Aromatic",
                                noun="Molecules",
                                verb="are",
                                label_type="binary",
                                calculator=aromatic,
                                mol_format="rd"),

                    Fingerprint(property="Atoms",
                                noun="Molecules",
                                verb="contain",
                                label_type="multiclass",
                                calculator=num_atoms,
                                mol_format="rd"),

                    Fingerprint(property="Rings",
                                noun="Molecules",
                                verb="contain",
                                label_type="multiclass",
                                calculator=num_rings,
                                mol_format="rd"),

                    Fingerprint(property="Branches",
                                noun="Molecules",
                                verb="contain",
                                label_type="multiclass",
                                calculator=num_branches,
                                mol_format="smiles"),

                    Fingerprint(property="Conjugated",
                                noun="Molecules",
                                verb="are",
                                label_type="binary",
                                calculator=conjugated,
                                mol_format="rd")]


"""Fingerprint(property="Electronic Spatial Extent",
            verb="",
            noun="Molecules",
            label_type="continuous",
            calculator=Electronic_Spatial_Extent,
            mol_format="smiles",
            file=True)]"""
