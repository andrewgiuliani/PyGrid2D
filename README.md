# PyGrid

This is a short Python code that generates two dimensional high order embedded boundary grids.  The code is fast, with most operations fully vectorized in Numpy.  There are a number of domains already implemented from [1], that include an annulus domain as well as the domain used in Ringleb flow.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/annulus.png" alt="annulus" width="300" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/ringleb.png" alt="annulus" width="300" >
</p>
<p align="center"> <i>High order embedded boundary annulus (left) and Ringleb (right) domains.  Each curved boundary segment has q+1 boundary points (red dots), where q = 3.</i> <p align="center">
  
## Goal
The goal of this work is to provide an easy to install, high order cut cell grid generator for simple boundary geometries that have a functional representation. PyGrid assumes that each cut cell has only one or two, possibly curved, edges associated to the embedded boundary. 

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/annulus_zoom.png" alt="annulus" width="250" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/ringleb_zoom.png" alt="annulus" width="250" >
</p>
<p align="center"> <i>Zooms onto cut cells of the annulus and Ringleb meshes.  Cut cells can have either one (left) or two curved edges (right).  Cut cells with two curved edges can occur on sharp corners of the embedded boundary.</i> <p align="center">
  
We do not allow split or tunneled cells for ease of code development.  This mesh generator does not handle mesh degenerecies, nor can it handle all types of cut cells.  

## How does it work?
1. The user provides a function, `in_domain`, that maps a spatial coordinate (x,y) to 1 if the point lies inside the domain or 0 if it does not.  Using this function, the regular grid points that lie in and out of the domain are determined.  For example, this is done in the figure below on the Ringleb domain.  The regular grid points that lie in the domain are shown in orange, while the regular grid points that are outside the domain are shown in black.

2.  When a regular grid point that lies outside the domain is adjacent another that is inside the domain, this means that between these two points, the embedded boundary crosses a Cartesian grid line.  Using the method of bisection, the precise point of intersection is computed.  Right now, the code assumes that there is only one point of intersection.

3. After these intersection points are computed, the cut cells are assembled.

4. If the user requests curved edges, then additional vertices that are approximately uniformly spaced (in arclength) along the embedded boundary are computed.  For circular embedded boundaries, we have explicit formulae to accomplish this.  For the Ringleb domain, this is done using the method of bisection.

5.  The cut cells are then written to file in the `.ply` format.  See here for more information http://paulbourke.net/dataformats/ply/.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/gridgen.png" alt="annulus" width="600" > 
</p>
<p align="center"> <i> Regular grid points that lie inside and outside the domain are respectively orange and black, computed in step 1.  Irregular grid points on the embedded boundary are red, computed in step 2.</i> <p align="center">

## 🏗&nbsp; Grid generation

The different arguments that `PyGrid.py` accepts are explained below:


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

## Contact
For help running the code, or any other questions, send me an email at
`giuliani AT cims DOT nyu DOT edu`

## Citing
If you find this code useful in your work, you can cite it with

Giuliani, Andrew and Berger, Marsha. "A state redistribution method for discontinuous Galerkin methods on curvilinear embedded boundary grids".

## 📓&nbsp; License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## References
[1] Giuliani, Andrew. "A two-dimensional stabilized discontinuous Galerkin method on curvilinear embedded boundary grids".

