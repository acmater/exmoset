from .general import general_fingerprints
from .containsatom import atom_fingerprints
from .containsbond import bond_fingerprints
from .substructure import substructure_fingerprints
from .fingerprint import Fingerprint

fingerprints = general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints
