# Molsim

Package to automate the identification of molecular similarity given an arbitrary set
of molecules and associated functions to calculate the value of particular properties.

This method is guided by the underlying principle of finding labels with minimum entropy
across the provided set of molecules. Intuitively these labels represent a "pure" description
of the space of interest, with the label adopting almost exclusively a single value.

Extending this approach to continuous labels is a future task with the discussion about
how this could be done provided below.

## API
The Similarity_Analysis Class handles the analysis of the molecular subset with a provided set
of properties. The properties were built with extensability in mind. The API for a property is
defined by the Property abstract class with its associated methods.

Property instances are further decomposed into three separate subclasses - atom, bond, and molecule
properties. All of these utilize the API provided by property, and the information provided
by similarity analysis is sorted in accordance with that.

### TODO

Figure out how to properly set abstract method in properties class to require that all
subclasses have the .values property without creating and AttributeError.

I would like to extend this code to continuous descriptors such as electronic spatial extent or something similar.
I think however that doing that would require implementing https://en.wikipedia.org/wiki/Limiting_density_of_discrete_points.
Which would be very tricky.

Additionally, these continuous decisions need to be taken against some standard I believe. What exactly
this standard should be is difficult.

The reason for this is that as they are not discrete, they will never have multiple instances of the same class, so the notion
of "purity" that is intuitive for discrete distributions is lost in the continuous case.

There are some alternatives https://thomasberrett.github.io/EntropyEstimation.pdf and https://github.com/paulbrodersen/entropy_estimators.

But the best formulation may be using KL-Div between the original distribution (in the entire dataset), and the distribution in the subset of interest.

1. Add the ability to print out molecules that could not be converted in Chem.mol objects
2. Another nice feature would be the capacity to select two distributions and identify what properties are different between them.
3. Add substructure properties which allow you to specify a functional group using SMARTS syntax and check to see if it is present in a particular set of molecules
   1. Is it possible to automate the identification of molecular features?
4. Speed it up. The code was written to get a working prototype as quickly as possible, but the current codebase is slow to execute and their should be some ways to dramatically accelerate it.
