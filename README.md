# openfoamautomation

## GOAL:
Use Genetic Algorithms to find the best scale, AoA and x y z position of multi step for a formula wing using CFDs

## Tech stack:
- Python - Language
- ANSA - Preprocessing
- Ansys Fluent (in replacement for ANSA temporary) - Preprocessing
- OpenFOAM - CFD 
- PyGAD - Genetic Library

## How to get started:
1. Install openfoam within your WSL env using sudo
2. Then clone this repository into your WSL env
3. "sudo apt install git" in your WSL env
4. Set your username and email: git config --global user.name "username" + git config --global user.email "email@blah.ca"
5. Make an ssh key using your email: ssh-keygen -t ed25519 -C "email@blah.ca"
6. Don't need a passphrase if you don't want to
7. Find your SSH using: cat ~/.ssh/id_ed25519.pub
8. Add the SSH key into your github in settings
9. Clone it into your WSL using: git clone git@github.com:username/repostiroyname.git
10. GO crazy

## Workflow (Draft)

- PyGAD puts in boundary conditions (Scale, AoA, x y z position of airfoils) then
- ANSA Pre Processing (Geometry cleanup, zone definition, meshing) then
- OpenFOAM (grabs .msh file) and uses their solver then 
- PyGAD grabs results from OpenFOAM then repeat step 1
- PyGAD spits out the best setup 

## Notes:

- Using Ansys instead of ANSA for now, until I get access + documentation for ANSA

