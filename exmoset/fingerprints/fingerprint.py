class Fingerprint():
    """
    Fingerprint that describes a particular property in accordance with a set of properties

    Parameters
    ----------
    property : str
        Name of the property.

    noun : "atom","bond","molecule","substructure"
        The noun in which the property should be calculated.

    label_type : "binary","multiclass","continuous"
        The type of label distribution that is expected - determines how the
        property is described and how the entropy is calculated.

    calculator : <function>
        A python function that will be called during the class determined by the noun
        to determine the value of the property.

    mol_format : str
        Which molecular format will be accessed from the Molecule object.
    """
    def __init__(self,property,
                      verb,
                      noun,
                      label_type,
                      calculator,
                      mol_format,
                      sensitivity=0.1):
        assert label_type in ["binary", "multiclass","continuous"], "Not a valid label type."
        self.property    = property
        self.verb        = verb
        self.noun        = noun
        self.label_type  = label_type
        self.calculator  = calculator
        self.mol_format  = mol_format
        self.sensitivity = sensitivity

    def __repr__(self):
        return f"Fingerprint({self.property},{self.verb},{self.noun},{self.label_type},{self.calculator},{self.mol_format},{self.sensitivity})"

    def __str__(self):
        return f"Fingerprint\n\t" + "\n\t".join([f"Property : {self.property:^11}",
                                                 f"Verb : {self.verb:^24}",
                                                 f"Noun: {self.noun:^18}",
                                                 f"Label Type: {self.label_type:^11}",
                                                 f"Calculator : {self.calculator}",
                                                 f"Mol Format : {self.mol_format:^8}"])

    def summary(self,val,entropy,sensitivity=None,unimportant_label=False):
        """
        A binary label summary.

        Grammatical Structure
        <noun> <verb> <property>

        The description is generated in accordance with the different language components.
        If the value is zero, then the statement is negated using a method in the label superclass.
        """
        print(val,entropy,sensitivity)
        if not sensitivity:
            sensitivity = self.sensitivity
        if entropy < sensitivity:

            if self.label_type == "binary":
                description = f"{self.noun} {self.verb} {self.property}"
                return description if val > 0.5 else self.neg(description)
            elif self.label_type == "multiclass":
                description = f"{self.noun} {self.verb} {val} {self.property}"
                return description
            else:
                return f"{self.noun} share a similar {self.property} centred around {val:.4f}"
        elif unimportant_label:
            return f"{self.property} not meaningful"

    @staticmethod
    def neg(string):
        if "are" in string:
            return string.replace("are", "are not")
        elif "contain" in string:
            return string.replace("contain", "do not contain")

if __name__ == "__main__":
    a = Fingerprint(property="Contains C",
                    noun="molecule",
                    label_type="binary",
                    calculator="add",
                    mol_format="smiles")
    print(a)
