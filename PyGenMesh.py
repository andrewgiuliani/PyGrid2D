import numpy as np
import ipdb
import sys
import argparse
import fbody as fb
import bisection as bs
import plot_mesh as pm

import ipdb



parser = argparse.ArgumentParser()
parser.add_argument("-Nx", "--Nx", type=int, help="number of elements in x-direction")
parser.add_argument("-Ny", "--Ny", type=int, help="number of elements in y-direction")
parser.add_argument("-left", "--LEFT", type=float, help="left endpoint", default = 0)
parser.add_argument("-right", "--RIGHT", type=float, help="right endpoint", default = 0)
parser.add_argument("-bottom", "--BOTTOM", type=float, help="h endpoint", default = 0)
parser.add_argument("-top", "--TOP", type=float, help="top endpoint", default = 0)
parser.add_argument("-fbody","--fbody", type=int, help="geometry type")
args = parser.parse_args()

Nx = args.Nx
Ny = args.Ny
bid = args.fbody
dom = np.array([args.LEFT, args.BOTTOM, args.RIGHT, args.TOP])

X = np.linspace(dom[0], dom[2], Nx+1)
Y = np.linspace(dom[1], dom[3], Ny+1)
XX,YY = np.meshgrid(X,Y)
XX = np.transpose(XX)
YY = np.transpose(YY)


vert_in = fb.fbody(XX,YY,bid)
idx_vert_in = np.where(vert_in)
num = vert_in[0:-1,0:-1] + vert_in[0:-1,1:] + vert_in[1:,0:-1] + vert_in[1:,1:]
irr = np.equal(num,4).astype(int)


cuts = np.logical_and(np.greater(num , 0) , np.less(num, 4) ).astype(int)
num_cuts = np.sum(cuts)

h = np.where( np.logical_xor(vert_in[0:-1,:],vert_in[1:,:]) ) 
v   = np.where( np.logical_xor(vert_in[:,0:-1],vert_in[:,1:]) ) 

xh1 = XX[h[0]  , h[1]]
xh2 = XX[h[0]+1, h[1]]
xh  = bs.bisection(xh1, xh2, lambda xin : fb.fbody(xin,YY[h],bid) )  

yv1 = YY[v]
yv2 = YY[v[0], v[1]+1]
yv  = bs.bisection(yv1, yv2, lambda yin : fb.fbody(XX[v],yin,bid) )  

whole_idx = np.where(irr)
cut_idx   = np.where(cuts)

has_bot    = np.logical_and(vert_in[cut_idx[0]   , cut_idx[1]  ],vert_in[cut_idx[0]+1 , cut_idx[1]  ])
has_right  = np.logical_and(vert_in[cut_idx[0]+1 , cut_idx[1]  ],vert_in[cut_idx[0]+1 , cut_idx[1]+1])
has_top    = np.logical_and(vert_in[cut_idx[0]   , cut_idx[1]+1],vert_in[cut_idx[0]+1 , cut_idx[1]+1])
has_left   = np.logical_and(vert_in[cut_idx[0]   , cut_idx[1]+1],vert_in[cut_idx[0]   , cut_idx[1]  ])

# make grid of irregular vertex indices
bot_idx = np.cumsum(has_bot)-1
right_idx = np.cumsum(has_right)-1 + bot_idx[
top_idx = np.cumsum(has_top)-1
left_idx = np.cumsum(has_left)-1


num_vertices = np.sum(vert_in) + xh.size + yv.size

# assemble all the vertices in the grid
v1 = whole_idx
v2 = (v1[0]+1, v1[1]  )
v3 = (v1[0]+1, v1[1]+1)
v4 = (v1[0]  , v1[1]+1)

# grid of indices assigned to Cartesian vertices
vertex_idx = np.cumsum(np.ravel(vert_in)).reshape(Nx+1,Ny+1)-1

# compress out the unused Cartesian vertices
regular_vertices = np.hstack( (
    np.compress(np.ravel(vert_in).astype(bool),np.ravel(XX))[:,None], 
    np.compress(np.ravel(vert_in).astype(bool),np.ravel(YY))[:,None]
    ) )
vertices = np.vstack( ( 
                     regular_vertices,
                     np.hstack((xh[:,None], YY[h][:,None] )),
                     np.hstack((XX[v][:,None], yv[:,None]))
                     ) )









# assemble the whole cells
cells_whole = np.hstack( (
                          vertex_idx[v1][:,None],
                          vertex_idx[v2][:,None],
                          vertex_idx[v3][:,None],
                          vertex_idx[v4][:,None]
                         ) )



# assemble the cut cells
v1 = cut_idx
v2 = (v1[0]+1, v1[1]  )
v3 = (v1[0]+1, v1[1]+1)
v4 = (v1[0]  , v1[1]+1)

v1_idx  = np.where( vert_in[v1][:,None],  vertex_idx[v1][:,None], -np.ones((num_cuts,1)).astype(int)  )
v12_idx = np.where(has_bot,bot_idx,-np.ones(num_cuts) )[:,None]
v2_idx  = np.where( vert_in[v2][:,None],  vertex_idx[v2][:,None], -np.ones((num_cuts,1)).astype(int)  )
v23_idx = np.where(has_right,right_idx,-np.ones(num_cuts) )[:,None]
v3_idx  = np.where( vert_in[v3][:,None],  vertex_idx[v3][:,None], -np.ones((num_cuts,1)).astype(int)  )
v34_idx = np.where(has_top,top_idx,-np.ones(num_cuts) )[:,None]
v4_idx  = np.where( vert_in[v4][:,None],  vertex_idx[v4][:,None], -np.ones((num_cuts,1)).astype(int)  )
v41_idx = np.where(has_left,left_idx,-np.ones(num_cuts) )[:,None]

cut_vertices = np.hstack( ( v1_idx,
                            v12_idx,
                            v2_idx,
                            v23_idx,
                            v3_idx,
                            v34_idx,
                            v4_idx,
                            v41_idx)
                        )
cut_nv = np.sum( cut_vertices != -1, axis = 1)



ipdb.set_trace(context=21)

# reorder vertices counterclockwise
#pm.plot_mesh(vertices,cells_whole)



#import matplotlib.pyplot as plt
#plt.scatter(vertices[:,0], vertices[:,1])
#plt.grid()
#plt.show()

