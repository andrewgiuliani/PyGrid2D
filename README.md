# PyGrid

This is a short Python code that generates two dimensional high order embedded boundary grids.  The code is fast, with most operations are fully vectorized in Numpy.  There are a number of domains already implemented from [1], that include an annulus domain as well as the domain used in Ringleb flow.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/annulus.png" alt="annulus" width="300" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/ringleb.png" alt="annulus" width="300" >
</p>
<p align="center"> <i>High order embedded boundary annulus (left) and Ringleb (right) domains.  Each curved boundary segment has q+1 boundary points (red dots), where q = 3.</i> <p align="center">
  
## Goal
The goal of this work is to provide an easy to install, high order cut cell grid generator for simple boundary geometries that have a functional representation. PyGrid assumes that each cut cell is assumed to have only one or two, possibly curved, edges associated to the embedded boundary. 

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/annulus_zoom.png" alt="annulus" width="250" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/ringleb_zoom.png" alt="annulus" width="250" >
</p>
<p align="center"> <i>Zooms onto cut cells of the annulus and Ringleb meshes.  Cut cells can have either one (left) or two curved edges (right).  Cut cells with two curved edges can occur on sharp corners of the embedded boundary.</i> <p align="center">
  
We do not allow split or tunneled cells for ease of code development.  This mesh generator does not handle mesh generecies, nor can it handle all types of cut cells.  

## 🏗&nbsp; Grid generation

The different arguments that `gengrid.py` accepts are explained below:


-Nx

number of cells on the grid in the x-direction

-Ny

number of cells on the grid in the y-direction

-fbody

the ID of the embedded boundary to be used.
* `0`: quarter annulus, for the supersonic vortex problem,
* `1`: annulus,
* `2`: Ringleb,
* `3`: random blobs.

-plot

display the generated cut cell grid

-q

add curved edges on the cut cells with q+1 interpolation points with the embedded boundary

For example, the Ringleb domain can be generated and plotted by calling
```
python PyGrid.py -Nx 10 -Ny 10 -fbody 2 -q 3 -plot
```
## Dependencies
For fancy command line output
```
pip install rich
```



## 📓&nbsp; License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## References
[1] Giuliani, Andrew and Berger, Marsha. "A state redistribution method for discontinuous Galerkin methods on curvilinear embedded boundary grids".

