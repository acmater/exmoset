from .general import general_fingerprints
from .containsatom import atom_fingerprints
from .containsbond import bond_fingerprints
from .substructure import substructure_fingerprints
from .fingerprint import Fingerprint
from .template_fingerprint import gen_fp, template_fp, mof_template

# This is where the default fingerprints are defined
fingerprints = general_fingerprints + atom_fingerprints + bond_fingerprints + substructure_fingerprints
