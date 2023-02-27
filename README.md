# Ch393 IS Force-Field Optimizer 
The project takes in an input file in .xyz and returns the optimized structure in .xyz with the associated energy value.

## Overall Strategy
1.	determine the xyz of each Ar atom by parsing the XYZ file.
2.	determine the distance (r) between the two atoms as a function of x1, y1, z1, x2, y2, y3.
3.	determine the potential using the distance in 2.
4.	update the XYZ coordinates of the two atoms as guided by the Nelder Mead algorithm.

## Progrgramming Tasks:
- [ ] Create a XYZ file of 1 Ar atom
- [ ] Create a XYZ file of 2 Ar atoms
- [ ] Parse the XYZ file
- [ ] Determine the distance between the 2 Ar atoms
- [ ] Determine LP potential energy between the 2 Ar atoms
- [x] Implement Lennard-Jones potential for Argon
- [x] Plot the potential energy curve
- [x] Determine the gradient using 

## Tools
Visual Studio Shortcuts:
https://code.visualstudio.com/shortcuts/keyboard-shortcuts-macos.pdf

## Keyboard Shortcuts:
`CTRL` + ```

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
