# openfoamautomation

## Usage:
Automated Genetic Algorithm to find optimal AoA and scale for airfoils using OpenFOAM CFD

## GOAL:
Use Genetic Algorithms to find the best scale, AoA and x y z position of multi-step for a formula wing using CFDs

## Tech stack:
- Python - Language

- OpenFOAM - CFD 
- PyGAD - Genetic Library
- ANSA - Preprocessing (Future Implementation)

## Workflow

2D Airfoil coords (.txt) -> Initial .env with BCs, including domain size + cells, flow conditions, turbulence models, air density, parallel core count, etc -> Run + Setup OpenFOAM CFDs w/ PyGAD

## Notes:

- Using OpenFOAM SnappyHexMesh

