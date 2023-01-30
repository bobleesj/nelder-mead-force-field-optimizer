# Force-Field Optimizer - Branch:
The project takes in an input file in .xyz and returns the optimized structure in .xyz with the associated energy value.

Branch is updated once again.



## To run

```
docker build -t main-python .        
docker run --name main-python-app main-python
```

```
docker rm main-python-app  # remove the container
docker build -t main-python . # re-build with the updated Python file
docker run --name main-python-app main-python # run
```