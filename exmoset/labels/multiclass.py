from exmoset.abstract import Label

class Multiclass(Label):
    """
    Class that handles the manipulation of multiclass molecular labels.
    """
    def __init__(self,name,
                      values,
                      context='atom'):
        self.name    = name
        self.values  = values
        self.av      = np.round(np.mean(values))

    def summary(self,verbose=False):
        if self.context == "atom":
            return f"Atoms share a similar {self.name} centred around {self.av}"
        elif self.context == "bonds":
            return f"Bonds share a similar {self.name} centred around {self.av}"
        elif self.context == "molecule":
            return f"Molecules share a similar {self.name} centred around {self.av}"
        else:
            return f"{self.context} {self.name} at {self.av}"
