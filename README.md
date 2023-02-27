# Ch393 IS Force-Field Optimizer 
The project takes in an input file in .xyz and returns the optimized structure in .xyz with the associated energy value.

## Overall Strategy
1.	determine the xyz of each Ar atom by parsing the XYZ file.
2.	determine the distance (r) between the two atoms as a function of x1, y1, z1, x2, y2, y3.
3.	determine the potential using the distance in 2.
4.	update the XYZ coordinates of the two atoms as guided by the Nelder Mead  algorithm.

## Further Features
- [ ] Optimize more than 2 Ar atoms

## Progrgramming Tasks:
- [ ] Randomize the initial positions of atoms
- [x] Refactor the code
- [x] Save the optimized coordinates into a xyz file.
- [x] Apply Nelder–Mead method to determine optimal coordinates
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

## Keyboard Shortcuts:
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