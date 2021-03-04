# EXplainable MOlecular SETs

### Considering possible namechange to the one above.

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties.

This method is guided by the underlying principle of finding labels with minimum entropy
across the provided set of molecules. Intuitively these labels represent a "pure" description
of the space of interest, with the label adopting almost exclusively a single value.

## API
The Similarity_Analysis Class handles the analysis of the molecular subset with a provided set
of properties. The properties were built with extensability in mind. The API for a property is
defined by the Property abstract class with its associated methods.

The majority of the work is handled by the property metaclass, which contains the key entropy functionality
for both continuous and discrete cases. This is called by all of its subclasses to determine the importance
of each of the properties that they are assessing.

Property instances are further decomposed into three separate subclasses - atom, bond, and molecule
properties. All of these utilize the API provided by property, and the information provided
by similarity analysis is sorted in accordance with that.

### TODO

I would like to extend this code to continuous descriptors such as electronic spatial extent or something similar.
I think however that doing that would require implementing https://en.wikipedia.org/wiki/Limiting_density_of_discrete_points.
Which would be very tricky.

Additionally, these continuous decisions need to be taken against some standard I believe. What exactly
this standard should be is difficult.

The reason for this is that as they are not discrete, they will never have multiple instances of the same class, so the notion
of "purity" that is intuitive for discrete distributions is lost in the continuous case.

There are some alternatives https://thomasberrett.github.io/EntropyEstimation.pdf and https://github.com/paulbrodersen/entropy_estimators.

But the best formulation may be using KL-Div between the original distribution (in the entire dataset), and the distribution in the subset of interest.

So the final major code refactorization that I want to consider now is a further decomposition of the problem.
Currently molprop, substructure, atom, and bond are the four main classes, but their functionality is restricted. Having these as their own superclasses that inherit from property and then spawning particular instances like atom_types would allow me to extend it further to other atomic properties such as charge, radical, hybridisation, etc.

2. Another nice feature would be the capacity to select two distributions and identify what properties are different between them.
3. Speed it up. The code was written to get a working prototype as quickly as possible, but the current codebase is slow to execute and their should be some ways to dramatically accelerate it.
6. Add code snippet examples for how to run

#### Notes
There are some portions of the code that look like they could be offloaded to the abstract class. In particular
```python
def summative_label
```
appears to duplicated. The reason for its repeated use is that each one is slightly different in the print statements that it provides. There may be a way to automate this by providing template strings, but currently that has not been examined.
