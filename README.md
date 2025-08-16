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

## How to get started:
1. Install openfoam within your Linux/WSL env using sudo
2. Then clone this repository into your env
3. "sudo apt install git" in your env
4. Set your username and email: git config --global user.name "username" + git config --global user.email "email@blah.ca"
5. Make an ssh key using your email: ssh-keygen -t ed25519 -C "email@blah.ca"
6. Don't need a passphrase if you don't want to
7. Find your SSH using: cat ~/.ssh/id_ed25519.pub
8. Add the SSH key into your github in settings
9. Clone it into your WSL using: git clone git@github.com:username/repostiroyname.git
10. GO crazy

## Notes:

- Using OpenFOAM SnappyHexMesh instead of ANSA for now, until I get access + documentation for ANSA

