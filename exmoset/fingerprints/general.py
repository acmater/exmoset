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
                                mol_format="rd"),

                    Fingerprint(name="Dipole Moment",
                                context="Molecules",
                                label_type="continuous",
                                calculator=Dipole_Moment,
                                mol_format="smiles",
                                file=True),

                    Fingerprint(name="Electronic Spatial Extent",
                                context="Molecules",
                                label_type="continuous",
                                calculator=Electronic_Spatial_Extent,
                                mol_format="smiles",
                                file=True)]
