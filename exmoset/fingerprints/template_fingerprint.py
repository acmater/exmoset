from .fingerprint import Fingerprint

template_fp = {"verb" : None,
               "noun" :"Molecules",
               "label_type" : "continuous",
               "calculator" : None,
               "mol_format" : "smiles"}

def gen_fp(prop):
    return Fingerprint(property=prop, **template_fp)
