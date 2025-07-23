# openfoamautomation

## GOAL: use generative algorithms to find the best scale, AoA and x y z position of multi step for a formula wing 
## Initialization: # of generations you want

Using Ansys instead of ANSA for now

PyGAD puts in boundary conditions (Scale, AoA, x y z position of airfoils) then
ANSA Pre Processing (Geometry cleanup, zone definition, meshing) then
OpenFOAM (grabs .msh file) and uses their solver then 
PyGAD grabs results from OpenFOAM then repeat step 1
PyGAD spits out the best setup 


## Tech stack:
ANSA Python Script Library
OpenFOAM Automation Script Library
PyGAD

gyatt
