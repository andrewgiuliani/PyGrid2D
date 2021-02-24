# PyGrid-2D


![build-and-test-status](https://github.com/andrewgiuliani/PyGrid-2D/workflows/Build%20&%20Test/badge.svg)



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
  
We do not allow split or tunneled cells for ease of code development.  This mesh generator does not handle mesh degeneracies, nor can it handle all types of cut cells.  

## How does it work?
1. The user provides a function, `in_domain`, that maps a spatial coordinate (x,y) to 1 if the point lies inside the domain or 0 if it does not.  Using this function, the regular grid points that lie in and out of the domain are determined.  For example, this is done in the figure below on the Ringleb domain.  The regular grid points that lie in the domain are shown in orange, while the regular grid points that are outside the domain are shown in black.

2.  When a regular grid point that lies outside the domain is adjacent another that is inside the domain, this means that between these two points, the embedded boundary crosses a Cartesian grid line.  Using the method of bisection, the precise point of intersection is computed.  On the Ringleb domain below, these points are plotted in red.  Right now, the code assumes that the boundary does not cross the grid line more than once.

3. After these intersection points are computed, the cut cells are assembled.

4. If the user requests curved edges, then additional vertices that are approximately uniformly spaced (in arclength) along the embedded boundary are computed.  For circular embedded boundaries, we have explicit formulae to accomplish this.  Again, this is done using the method of bisection.

5.  The cut cells are then written to file in the `.ply` format.  See here for more information http://paulbourke.net/dataformats/ply/.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/gridgen.png" alt="annulus" width="600" > 
</p>
<p align="center"> <i> Regular grid points that lie inside and outside the domain are respectively orange and black, computed in step 1.  Irregular grid points on the embedded boundary are red, computed in step 2.</i> <p align="center">
  
## 🌲🌳&nbsp; Into the weeds on curved edges
If curved boundaries are requested, the code computes interpolation points that are approximately equispaced in arclength using the method of bisection.  This approach relies on a one-dimensional parametrization of the boundary  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Gamma(s)&space;:&space;s&space;\rightarrow&space;(x(s),y(s))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Gamma(s)&space;:&space;s&space;\rightarrow&space;(x(s),y(s))" title="\Gamma(s) : s \rightarrow (x(s),y(s))" /></a>

Say that a cut cell's curved edge defined by the set of points <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;s&space;\in&space;[s_1,s_2]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;s&space;\in&space;[s_1,s_2]" title="s \in [s_1,s_2]" /></a>.

Define the function

<a href="https://www.codecogs.com/eqnedit.php?latex=f(s)&space;=&space;\frac{||\Gamma(s)&space;-&space;\Gamma(s_1))||_2}{||\Gamma(s_2)-\Gamma(s_1)||_2}-w" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f(s)&space;=&space;\frac{||\Gamma(s)&space;-&space;\Gamma(s_1))||_2}{||\Gamma(s_2)-\Gamma(s_1)||_2}-w" title="f(s) = \frac{||\Gamma(s) - \Gamma(s_1))||_2}{||\Gamma(s_2)-\Gamma(s_1)||_2}-w" /></a>

for some <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;w&space;\in&space;[0,1]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;w&space;\in&space;[0,1]" title="w \in [0,1]" /></a>.  This function is zero when s maps to a point that is exactly <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;w&space;||&space;\Gamma(s_2)&space;-&space;\Gamma(s_1)||_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;w&space;||&space;\Gamma(s_2)&space;-&space;\Gamma(s_1)||_2" title="w || \Gamma(s_2) - \Gamma(s_1)||_2" /></a> units away from <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Gamma(s_1)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Gamma(s_1)" title="\Gamma(s_1)" /></a>.

## 🔪&nbsp; The sharp bits
This code only handles cut cells on which the boundary enters and leaves on different faces.  As a result, it cannot deal with tunneled and split cells for the moment.  If these scenarios are detected, the code returns an error.  
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/tunneled.png" alt="tunneled" width="250" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/split.png" alt="split" width="250" >
</p>
<p align="center"> <i> On the left, the annulus domain is removed from the background grid, creating tunneled cut cells.  On the right, the complement of the annulus is removed from the background grid, creating split cut cells. </i> <p align="center">
  
Additionally, the code does not robustly (or gracefully) handle mesh degeneracies.  For example, if the embedded geometry lies directly on a Cartesian grid line cell, the code will produce erroneous output and may, or may not output an error.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/aligned.png" alt="aligned" width="250" > 
</p>
<p align="center"> <i> A square domain (black lines) is removed from the background grid (blue lines), where the domain is perfectly aligned with the Cartesian grid lines. </i> <p align="center">
  
We do not support the above scenarios for code simplicity, however, they must be addressed when moving to three dimensions and complex engineering geometries.

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

Giuliani, Andrew. "A two-dimensional stabilized discontinuous Galerkin method on curvilinear embedded boundary grids". [https://arxiv.org/abs/2005.05734](https://arxiv.org/abs/2005.05734)

## 📓&nbsp; License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## References
[1] Giuliani, Andrew. "A two-dimensional stabilized discontinuous Galerkin method on curvilinear embedded boundary grids". [https://arxiv.org/abs/2005.05734](https://arxiv.org/abs/2005.05734)

