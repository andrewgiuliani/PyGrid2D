import numpy as np
import ssv
import rh
import ringleb as rb
import supersonic as ssv
import annulus as rh

bnames = {
            0 :'annulus', 
            1 :'supersonic',
            2 :'ringleb'
            }


def fbody(X,Y,bid):
    if bid == 0:
        active = rh.in_domain(X,Y)
    elif bid == 1:
        active = ssv.in_domain(X,Y)
    elif bid == 2:
        active = rb.in_domain(X,Y)
    else:
        active = np.ones(X.shape)

    return active.astype(int)

def bname(bid):
    return bnames[bid]


def is_corner(irreg_vertices, vertices,bid):
    X1 = vertices[irreg_vertices[:,0],0]
    Y1 = vertices[irreg_vertices[:,0],1]

    X2 = vertices[irreg_vertices[:,1],0]
    Y2 = vertices[irreg_vertices[:,1],1]

    if bid == 2:
        return rb.is_corner(X1,Y1, X2,Y2)
    else:
        return np.zeros(irreg_vertices.shape[0])

def compute_curved1(edges,vertices, q, bid):
    extra_vertices = np.zeros(edges.shape[0],q-1,2)
    wgt = np.linspace(0.,1.,qq+1)
    X1 = vertices[edges[:,0],0]
    Y1 = vertices[edges[:,0],1]
    X2 = vertices[edges[:,1],0]
    Y2 = vertices[edges[:,1],1]
    for qq in range(1,qq+1):
        extra_vertices[:,qq-1]=curved_body_ids[bid](wgt,X1,Y1,X2,Y2)
