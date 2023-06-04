# Choosing between models and boundary descriptions in wound healing applications

This code is used to generate the simulations corresponding to Chapter 4 of my PhD thesis. This code uses Chaste, details provided below. Code in this repository is modified from https://github.com/jmosborne/TissueBoundaries (code developed by Domenic Germano, James Osborne and myself for a paper submitted to the Bulletin of Mathematical Biology).

## Chaste
The code in this repository should run using the core version of Chaste, which can be found [here](https://chaste.cs.ox.ac.uk/trac/wiki), along with documentation on how to install and run simulations. Note that Chaste is only fully supported on Linux/Unix systems, Windows and Mac OS users are encouraged to use [Docker](https://github.com/Chaste/chaste-docker). 
 
To use the code in this repository, clone this repository into the `projects` subdirectory of the Chaste code. 

## Repository structure
- `src/` additional class files required to run the simulations using Chaste
	* `GhostNodeRemovalModifier.hpp` to remove ghost nodes when the area of the their associated Voronoi polygon becomes too small
	* `VertexBoundaryRefinementModifier.hpp` to add and remove nodes in the smooth and curved VD model
	* `VoidAreaModifier.hpp` to track the size of the wound
-`test/` test files for simulation setup
	* `TestInternalVoid.hpp` to generate simulations
-`figures/` MATLAB code to generate area over time curves (all other figures were generated using Paraview)
	* `InternalVoidAreaPlots.m` to plot wound area over time curves

## Running simulations using Chaste
To run the code (compiling using scons):
```
~: cd [Chaste Directory]
~/Chaste: scons b=GccOpt co=1 projects/WoundHealingModelSelection/test/TestInternalVoid.hpp &> build_output.log
~/Chaste: ./projects/WoundHealingModelSelection/build/optimised/TestInternalVoidRunner &> run_output.log
```
 