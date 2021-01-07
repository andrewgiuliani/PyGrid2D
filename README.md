# PyGrid

This is a short Python code that generates two dimensional high order embedded boundary grids.  The code is fast, with most operations are fully vectorized in Numpy.  There are a number of domains already implemented from [1], that include an annulus domain as well as the domain used in Ringleb flow.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/annulus.png" alt="annulus" width="300" > &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
  <img src="https://github.com/andrewgiuliani/PyGrid/blob/main/images/ringleb.png" alt="annulus" width="300" >
</p>
<p align="center"> <i>High order embedded boundary annulus (left) and Ringleb (right) domains.  Each curved boundary segment has q+1 boundary points (red dots), where q = 3.</i> <p align="center">
  
## Goal
The goal of this work is to provide an easy to install, high order cut cell grid generator for simple boundary geometries that have a functional representation.  This mesh generator does not handle mesh generecies, nor can it handle all types of cut cells.


## 🏗&nbsp; Grid generation


## 📓&nbsp; License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## References
[1] Giuliani, Andrew and Berger, Marsha. "A state redistribution method for discontinuous Galerkin methods on curvilinear embedded boundary grids".

