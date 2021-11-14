import pytest
import numpy as np
import pygrid2d as pg2d


def load_vertices_elem(domain, Nx, Ny, q):
    f = open("./data/" + domain.name+"_"+str(Nx)+"_"+str(Ny)+"_q" + str(q)+".ply", 'r')
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline().split(" ")
    Nv = int(line[-1])

    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline()
    
    vertices = np.zeros( (Nv, 2) )
    for v in range(Nv):
        line =  f.readline().split(" ")
        line = [float(l.strip('\n')) for l in line]
        line = np.array(line)
        vertices[v,:] = line 
    
    return vertices

@pytest.mark.parametrize("Nx,Ny,q,bid", [
    (10,10,1,0),
    (20,20,2,0),
    (30,30,3,0),
    (40,40,4,0),
    (50,50,5,0),
    (10,10,1,1),
    (20,20,2,1),
    (30,30,3,1),
    (40,40,4,1),
    (50,50,5,1),
    (10,10,1,2),
    (20,20,2,2),
    (30,30,3,2),
    (40,40,4,2)
])
def test_PyGrid(Nx, Ny, q, bid):
    vertices, elem, domain = pg2d.PyGrid2D(Nx, Ny, False, q, bid)
    vertices_loaded = load_vertices_elem(domain, Nx, Ny, q)

    err = np.linalg.norm(vertices - vertices_loaded)
    assert err < 1e-12
