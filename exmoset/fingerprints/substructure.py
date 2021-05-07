"""
Prewritten fingerprints for substructure matching based properties.
"""

import numpy as np
from .fingerprint import Fingerprint
from rdkit import Chem

# Location that substructures are drawn from - https://daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

substructures = ['*=O',
                 '[NX3;H2,H1,H0;!$(NC=[!#6])]',  # Must not be bonded to anything other than carbon
                 '[OX2H]([#6X4])', # Must be connected to sp3
                 '*O*',
                 '[OD2]([#6X4])[#6X4]', # Must also be connected to sp3
                 'C#N',
                 '[#6]([CX3](=O)[OX2H1])', # Must be connected to another carbon, removes the possibility of carbamates or similar groups.
                 '[$([CX3]=[OX1]),$([CX3+]-[OX1-])]',
                 '[NX3][CX3](=[OX1])[#6]']

def contains(substructure):
    match = Chem.MolFromSmarts(substructure)
    def sub_contains(mol):
        return int(mol.HasSubstructMatch(match))
    return sub_contains

substructure_fingerprints = []

for substructure in substructures:
    substructure_fingerprints.append(Fingerprint(property=f"{substructure}",
                verb="contain",
                noun="Molecules",
                label_type="binary",
                calculator=contains(substructure),
                mol_format="rd"))
