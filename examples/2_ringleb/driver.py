#!/usr/bin/env python3
import pygrid2d as pg2d

Nx = 25
Ny = 25
plot_flag = True
q = 5   # q = 5 means the mesh generator will compute 6 points on the embedded boundary, resulting in a degree 5 polynomial boundary interpolant
bid = 2 # bodyID = 2 means Ringleb domain

vertices, cell_list, domain, mesh_data, face_data = pg2d.PyGrid2D(Nx, Ny, plot_flag, q, bid)
pg2d.output_unstructured(vertices, mesh_data, face_data, domain, Nx, Ny, q)
pg2d.output_ply(vertices, cell_list, domain, Nx, Ny, q)
