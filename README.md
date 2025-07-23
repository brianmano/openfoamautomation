# openfoamautomation


Need to:
1. Install openfoam within your WSL env using sudo
2. Then clone this repository into your WSL env
    2.1. "sudo apt install git" in your WSL env
    2.2. Set your username and email: git config --global user.name "username" + git config --global user.email "email@blah.ca"
    2.3. Make an ssh key using your email: ssh-keygen -t ed25519 -C "email@blah.ca"
        2.3.1. Don't need a passphrase if you don't want to
    2.4. Find your SSH using: cat ~/.ssh/id_ed25519.pub
    2.5. Add the SSH key into your github in settings
    2.6. Clone it into your WSL using: git clone git@github.com:username/repostiroyname.git
    2.7. GO crazy

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
