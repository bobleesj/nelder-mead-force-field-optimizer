# Ch393 IS Force-Field Optimizer 
The project takes in an input file in .xyz and returns the optimized structure in .xyz with the associated Lennard-Jones potential valuve given r, ε, and α.

## Overall Strategy for 

## Programming Tasks:
- [ ] Save the optimized file
- [x] Apply Neader-Mead for N atoms
- [x] Determien the number of LJ interactions for N atoms
- [x] Propose an implementation strategy to optimize more than 2 Ar atoms
- [x] Randomize initial positions of atoms
- [x] Refactor code into modular functions
- [x] Save the optimized coordinates into a xyz file.
- [x] Apply Nelder–Mead method to find optimized structures
- [x] Express distance as a function of 2 points of (x, y, z)
- [x] Express the coordinates as a parameters
- [x] Determine LP potential energy between the 2 Ar 
- [x] Create a XYZ file of 2 Ar atoms using Avogadro
- [x] Parse the XYZ file using the xyz_py library
- [x] Determine the distance between the 2 Ar atomsatoms
- [x] Implement Lennard-Jones potential for Argon
- [x] Plot the potential energy curve

## Toolsgit
Visual Studio Shortcuts:
https://code.visualstudio.com/shortcuts/keyboard-shortcuts-macos.pdf

## Visual Studio Keyboard Shortcuts:
`CTRL` + ``` - Toggle Terminal
`CTRL` + g - Go to Line

## To run the file
```
docker build -t main-python .        
docker run --name main-python-app main-python
```

```
docker rm main-python-app  # remove the container
docker build -t main-python . # re-build with the updated Python file
docker run --name main-python-app main-python # run
```