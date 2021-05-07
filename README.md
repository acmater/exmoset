# EXplainable MOlecular SETs

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties (label fingerprints).

# Installation
The easiest way to install is using pip:

```bash
pip install exmoset
```

## API
The `MolSpace` Class handles the analysis of a given molecular set in accordance with the list of fingerprints provided. The molecules can be passed to Molspace in any format, with additional conversions specified by the `mol_converters` argument.

```python
analysis = MolSpace(molecules,
                    fingerprints = fingerprints,
                    file="data/QM9_Data.csv",
                    mol_converters={"rd" : Chem.MolFromSmiles, "smiles" : str},
                    index_col="SMILES")
```

### Fingerprints
Fingerprints are a standardized way for Molspace to calculate the properties for each molecule it is analysing. Its arguments determine the grammatical structure of the label that will be produced (`property, noun and verb`), and a function to calculate the property (`calculator`) along with what molecular format this function works on (`mol_format`). The grammatical structure of the resulting labels is a work in progress, and may lead to some poor results that require further processing.

```python
def contains_C(mol):
      return 1 if C in mol else 0

contains_carbon = Fingerprint(property="Contains C",
                  verb="contain",
                  noun="Molecule",
                  label_type="binary",
                  calculator=contains_C,
                  mol_format="smiles")
```

### Molecule Converters
The mol_converters argument provides the means to transform each molecule into alternate representations. The argument is a dictionary with the following structure {Identifier : Function_that_will_convert} that is expanded in the following way:

```python
formats = {key : mol_converters[key](mol) for key in mol_converters.keys()} # Assigns each identifier to its assocaited representation by
self.Molecules.append(Molecule(mol, **formats)) # Unpacks the new formats as kwargs into the Molecule object
```
An example is provided below
```python
mol_converters = {"rd" : Chem.MolFromSmiles, "smiles" : str} # Will convert molecules provided as smiles strings into Chem.rd objects from RDKit and maintain the SMILES in the dataset as strings.
```

## Label Types
### Binary
Binary labels indicate the presence of absence of a particular element, bond type, or molecular feature (such as aromaticity). Simplest to calculate and best behaved with respect to the entropy estimators. Uses a discrete entropy estimator.

### Multiclass
Discrete labels where the value can be any integer. Examples include number of rings, number of atoms, or number of each type of bond. Uses a discrete entropy estimator.

### Continuous
Continuous labels where the value can be any real number. Examples include electronic spatial extent, dipole moment, and free energy. Uses the continuous entropy estimator

## References
Almost all code in this work is original and based on the work of the following authors.

The exception is continuous entropy estimation which is provided by Paul Broderson's entropy estimators package (https://github.com/paulbrodersen/entropy_estimators).
