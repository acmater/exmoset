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
    return np.array([bool(mol.GetAromaticAtoms())],dtype=np.int)
def num_atoms(mol):
    return np.array([len(mol.GetAtoms())])
def num_rings(mol):
    return np.array([mol.GetRingInfo().NumRings()])


general_fingerprints = [Fingerprint(property="Aromatic",
                                context="Molecules",
                                verb="are",
                                label_type="binary",
                                calculator=aromatic,
                                mol_format="rd"),

                    Fingerprint(property="Atoms",
                                context="Molecules",
                                verb="contain",
                                label_type="multiclass",
                                calculator=num_atoms,
                                mol_format="rd"),

                    Fingerprint(property="Rings",
                                context="Molecules",
                                verb="contain",
                                label_type="multiclass",
                                calculator=num_rings,
                                mol_format="rd"),

                    Fingerprint(property="Dipole Moment",
                                verb="",
                                context="Molecules",
                                label_type="continuous",
                                calculator=Dipole_Moment,
                                mol_format="smiles",
                                file=True)]

"""Fingerprint(property="Electronic Spatial Extent",
            verb="",
            context="Molecules",
            label_type="continuous",
            calculator=Electronic_Spatial_Extent,
            mol_format="smiles",
            file=True)]"""
