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

# grid of booleans that say if horizontal or vertical faces have irregular intersections
# h_bool[0][0] is TRUE ==> cell 0 has a bottom irreg
# h_bool[0][1] is TRUE ==> cell 0 has a top irreg
h_bool  = np.logical_xor(vert_in[0:-1,:],vert_in[1:,:])
v_bool  = np.logical_xor(vert_in[:,0:-1],vert_in[:,1:])
h_idx   = np.where( h_bool ) 
v_idx   = np.where( v_bool ) 

xh1 = XX[h_idx[0]  , h_idx[1]]
xh2 = XX[h_idx[0]+1, h_idx[1]]
xh  = bs.bisection(xh1, xh2, lambda xin : fb.fbody(xin,YY[h_idx],bid) )  

yv1 = YY[v_idx]
yv2 = YY[v_idx[0], v_idx[1]+1]
yv  = bs.bisection(yv1, yv2, lambda yin : fb.fbody(XX[v_idx],yin,bid) )  

whole_idx = np.where(irr)
cut_idx   = np.where(cuts)


has_bot  = h_bool[ cut_idx[0]  , cut_idx[1]  ]
has_top  = h_bool[ cut_idx[0]  , cut_idx[1]+1]
has_left = v_bool[ cut_idx[0]  , cut_idx[1]  ]
has_right= v_bool[ cut_idx[0]+1, cut_idx[1]  ]


# make grid of irregular vertex indices
h_vert_idx   = np.cumsum(np.ravel(h_bool)).reshape(h_bool.shape)-1
v_vert_idx   = np.cumsum(np.ravel(v_bool)).reshape(v_bool.shape)-1 + (h_vert_idx[-1,-1]+1)


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
                     np.hstack((xh[:,None], YY[h_idx][:,None] )),
                     np.hstack((XX[v_idx][:,None], yv[:,None]))
                     ) )









# assemble the whole cells
cells_whole = np.hstack( (
                          vertex_idx[v1][:,None],
                          vertex_idx[v2][:,None],
                          vertex_idx[v3][:,None],
                          vertex_idx[v4][:,None]
                         ) )



# assemble the cut cells

# shift the irregular vertices to be at the end of the vertices list
h_vert_idx = h_vert_idx + np.sum(vert_in)
v_vert_idx = v_vert_idx + np.sum(vert_in) 

v1 = cut_idx
v2 = (v1[0]+1, v1[1]  )
v3 = (v1[0]+1, v1[1]+1)
v4 = (v1[0]  , v1[1]+1)

v1_idx  = np.where( vert_in[v1][:,None],  vertex_idx[v1][:,None], -np.ones((num_cuts,1)).astype(int)  )
v12_idx = np.where(has_bot,h_vert_idx[v1],-np.ones(num_cuts) )[:,None]

v2_idx  = np.where( vert_in[v2][:,None],  vertex_idx[v2][:,None], -np.ones((num_cuts,1)).astype(int)  )
v23_idx = np.where(has_right,v_vert_idx[v2],-np.ones(num_cuts) )[:,None]

v3_idx  = np.where( vert_in[v3][:,None],  vertex_idx[v3][:,None], -np.ones((num_cuts,1)).astype(int)  )
v34_idx = np.where(has_top,h_vert_idx[v4],-np.ones(num_cuts) )[:,None]

v4_idx  = np.where( vert_in[v4][:,None],  vertex_idx[v4][:,None], -np.ones((num_cuts,1)).astype(int)  )
v41_idx = np.where(has_left,v_vert_idx[v1],-np.ones(num_cuts) )[:,None]

cut_vertices = np.hstack( ( v1_idx,v12_idx,v2_idx,v23_idx,v3_idx,v34_idx,v4_idx,v41_idx)  ).astype(int)
cut_nv = np.sum( cut_vertices != -1, axis = 1)

s_idx = np.argsort(cut_nv)
cut_nv = cut_nv[s_idx]
cut_vertices = cut_vertices[s_idx,:]


cut_cells = [None] * 4 # here, we only support cut cells with 3,4,5, or 6 vertices
for nv in range(3,6):
    idx = np.where(cut_nv == nv)[0]
    cells_nv = np.ravel(cut_vertices[idx,:])
    
    pos_vals = cells_nv > -1
    cell_nv = np.compress(pos_vals, cells_nv).reshape((-1,nv))
    cut_cells[nv-3] = cell_nv

cell_list = [cells_whole] + cut_cells
vertex_count = [4,3,4,5]


# compute mesh stats
def PolyArea(x,y):
    return 0.5*np.abs(np.sum(x * np.roll(y,1,axis=1),axis = 1)-np.sum(y * np.roll(x,1,axis=1), axis = 1))
area = np.zeros((0,1))
for cells,nv in zip(cell_list, vertex_count):
    xcoord = vertices[np.ravel(cells),0].reshape( (cells.shape[0], nv) )
    ycoord = vertices[np.ravel(cells),1].reshape( (cells.shape[0], nv) )
    ar = PolyArea(xcoord,ycoord)
    area = np.vstack( (area,ar[:,None]) )


reg_vol = area[0]
min_vol_frac = np.min(area) / reg_vol
str_whole = str(cells_whole.shape[0])
str_cut = str(num_cuts)
str_tot = str(cells_whole.shape[0]+num_cuts)
str_min_vol = '%.16e' % min_vol_frac[0]

from rich.console import Console
from rich.table import Column, Table
console = Console()

table = Table(show_header = False)
table.add_column("Name", style="dim", width=15)
table.add_column("Value",justify="right")
table.add_row(
    "whole cells", str_whole
)
table.add_row(
    "cut cells", str_cut 
        )
table.add_row(
    "total cells", str_tot 
)
table.add_row(
    "min. vol. frac.", str_min_vol 
)
console.print(table)

#str_whole = str(cells_whole.shape[0])
#str_cut = str(num_cuts)
#str_tot = str(cells_whole.shape[0]+num_cuts)
#stats = 'whole cells \t {whole} \ncut cells \t {cut} \ntotal number \t {tot} '.format(whole = str_whole, cut = str_cut, tot = str_tot) 
#print('**** MESH METADATA ****')
#print(stats) 
#print('***********************')
# reorder vertices counterclockwise
#pm.plot_mesh(vertices,cell_list,vertex_count)


#ipdb.set_trace(context=21)

#import matplotlib.pyplot as plt
#plt.scatter(vertices[:,0], vertices[:,1])
#plt.grid()
#plt.show()

