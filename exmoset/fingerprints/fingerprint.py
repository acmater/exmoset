class Fingerprint():
    """
    Fingerprint that describes a particular property in accordance with a set of properties

    Parameters
    ----------
    name : str
        Name of the property.

    context : "atom","bond","molecule","substructure"
        The context in which the property should be calculated.

    label_type : "binary","multiclass","continuous"
        The type of label distribution that is expected - determines how the
        property is described and how the entropy is calculated.

    calculator : <function>
        A python function that will be called during the class determined by the context
        to determine the value of the property.

    mol_format : str
        Which molecular format will be accessed from the Molecule object.

    file : str, default=None
        A file that can be optionally indexed by the fingerprint method.
    """
    def __init__(self,name,
                      context,
                      label_type,
                      calculator,
                      mol_format,
                      file=None):
        self.name       = name
        self.context    = context
        self.label_type = label_type
        self.calculator = calculator
        self.mol_format = mol_format
        self.file       = file

    def __repr__(self):
        return f"Fingerprint({self.name},{self.context},{self.label_type},{self.calculator},{self.mol_format})"

    def __str__(self):
        return f"Fingerprint\n\t" + "\n\t".join([f"Name : {self.name:^24}",
                                                 f"Context: {self.context:^18}",
                                                 f"Label Type: {self.label_type:^11}",
                                                 f"Calculator : {self.calculator:^5}",
                                                 f"Mol Format : {self.mol_format:^8}"])

if __name__ == "__main__":
    a = Fingerprint(name="Contains C",
                    context="molecule",
                    label_type="binary",
                    calculator="add",
                    mol_format="smiles")
    print(a)
