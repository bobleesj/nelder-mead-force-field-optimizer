# Ch393 IS Force-Field Optimizer 
The project takes in an input file in .xyz and returns the optimized structure in .xyz with the associated Lennard-Jones potential valuve given r, ε, and α.

## Overall Strategy for 3D geometry optimization of 2 Argon atoms (Implemented)
1.	determine the xyz coordinates of each Ar atom by parsing the XYZ file.
2.	determine the distance (r) between the two atoms as a function of x1, y1, z1, x2, y2, y3.
3.	determine the potential using the distance in 2.
4.	update the XYZ coordinates of the two atoms as guided by the Nelder-Mead algorithm.

## Overall Strategy for 3D geometry optimization of 3 Argon atoms (To be implemented)
1. determine the xyz coordinates of each Ar atoms by parsing the XYZ file.
2. determine the distances (r1, r2, r3) of Ar_1-Ar_2, Ar-2-Ar-3, and Ar3-Ar-1.
3. determine the potential as a function of r1, r2, r3. Note each r is a function of x, y, and z.
4. update the 9 coordinates towards an optimal using the Nelder-Mead algorithm.

## Overall Strategy for 3D geometry optimization of N Argon atoms (To be implemented)
TBD

## Overall Strategy for 

## Progrgramming Tasks:
- [x] Propose an implementation strategy to optimize more than 2 Ar atoms
- [x] Randomize initial positions of atoms
- [x] Refactor code into modular functions
- [x] Save the optimized coordinates into a xyz file.
- [x] Apply Nelder–Mead method to determine 3D coordinates
- [x] Express distance as a function of 2 points of x, y, z
- [x] Express the coordinates as a parameters
- [x] Determine LP potential energy between the 2 Ar 
- [x] Create a XYZ file of 2 Ar atoms using Avagadro
- [x] Parse the XYZ file using the xyz_py library
- [x] Determine the distance between the 2 Ar atomsatoms
- [x] Implement Lennard-Jones potential for Argon
- [x] Plot the potential energy curve
- [x] Determine the gradient using 

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