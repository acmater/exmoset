from .fingerprint import Fingerprint

template_fp = {"verb" : None,
               "noun" :"Molecules",
               "label_type" : "continuous",
               "calculator" : None,
               "mol_format" : "smiles"}

mof_template = {"verb" : None,
               "noun" :"MOFs",
               "label_type" : "continuous",
               "calculator" : None,
               "mol_format" : "mofid"}


def gen_fp(prop,template=template_fp):
    return Fingerprint(property=prop, **template)
