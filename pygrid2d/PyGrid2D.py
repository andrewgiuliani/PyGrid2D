import numpy as np
import sys
import argparse
from .bisection import *
from .plot_mesh import *
from .output import *
from .domain import Domain
from .ringleb import Ringleb
from .annulus import Annulus
from .annulus_acoustics import Annulus_acoustics
from .supersonic import Supersonic
from .blobs import Blobs
from .circles import Circles
from .channel import Channel
from .doublemach import DoubleMach
from rich.console import Console
from rich.table import Column, Table

def PyGrid2D(Nx, Ny, plot_flag, q, bid):

    console = Console()
    
    if bid == 0:
        domain = Annulus()
    elif bid == 1:
        domain = Supersonic()
    elif bid == 2:
        domain = Ringleb()
    elif bid == 3:
        domain = Blobs()
    elif bid == 4:
        domain = Annulus_acoustics()
    elif bid == 5:
        domain = Circles()
    elif bid == 6:
        domain = Channel()
    elif bid == 7:
        domain = DoubleMach()
    else:
        quit()
    
    
    dom = np.array([domain.left, domain.bottom, domain.right, domain.top])
    console.print("Embedded grid input", style="bold red")
    table = Table(show_header = False)
    table.add_column("Name", style="dim", width=15)
    table.add_column("Value",justify="right")
    table.add_row("(Nx,Ny)", '({},{})'.format(Nx,Ny))
    table.add_row("Domain", '[{},{}] x [{},{}]'.format(dom[0],dom[2], dom[1], dom[3] ))
    table.add_row("Grid geometry", domain.name )
    console.print(table)
    
    
    
    console.print("Generating the grid...", style="bold blue")
    X = np.linspace(dom[0], dom[2], Nx+1)
    Y = np.linspace(dom[1], dom[3], Ny+1)
    XX,YY = np.meshgrid(X,Y)
    XX = np.transpose(XX)
    YY = np.transpose(YY)
    
    
    vert_in = domain.in_domain(XX,YY)
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
    
    console.print("Computing the x-intersections...", style="bold blue")
    xh1 = XX[h_idx[0]  , h_idx[1]]
    xh2 = XX[h_idx[0]+1, h_idx[1]]
    xh  = bisection(xh1, xh2, lambda xin : (domain.in_domain(xin,YY[h_idx])-0.5) )  
    
    console.print("Computing the y-intersections...", style="bold blue")
    yv1 = YY[v_idx]
    yv2 = YY[v_idx[0], v_idx[1]+1]
    yv  = bisection(yv1, yv2, lambda yin : (domain.in_domain(XX[v_idx],yin)-0.5) )  
    
    whole_idx = np.where(irr)
    cut_idx   = np.where(cuts)
    
    
    has_bot  = h_bool[ cut_idx[0]  , cut_idx[1]  ]
    has_top  = h_bool[ cut_idx[0]  , cut_idx[1]+1]
    has_left = v_bool[ cut_idx[0]  , cut_idx[1]  ]
    has_right= v_bool[ cut_idx[0]+1, cut_idx[1]  ]
    
    
    # make grid of irregular vertex indices
    h_vert_idx   = np.cumsum(np.ravel(h_bool)).reshape(h_bool.shape)-1
    v_vert_idx   = np.cumsum(np.ravel(v_bool)).reshape(v_bool.shape)-1 + (h_vert_idx[-1,-1]+1)
    
    num_regular_vertices = np.sum(vert_in)
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
    
    
    
    
    
    console.print("Shifting the cut cell vertices...", style="bold blue")
    
    s_idx = np.argsort(cut_nv)
    cut_nv = cut_nv[s_idx]
    cut_vertices = cut_vertices[s_idx,:]
    cut_idx = (cut_idx[0][s_idx], cut_idx[1][s_idx])
    
    irreg_edges = np.zeros( ( num_cuts, 2 ) ).astype(int)
    
    
    
    min_nv = np.min(cut_nv)
    max_nv = np.max(cut_nv)
    cut_cells = [np.array([]) for vv in np.arange(min_nv,max_nv+1)]
    cell_list_idx = [ () for vv in np.arange(min_nv,max_nv+1)]
    
    count = 0
    for nv in np.arange(min_nv,max_nv+1):
        idx = np.where(cut_nv == nv)[0]
        cells_nv = np.ravel(cut_vertices[idx,:])
        
        pos_vals = cells_nv > -1
        cell_nv = np.compress(pos_vals, cells_nv).reshape((-1,nv))
        
        # this is another shift to make sure first and last vertices are edge vertices
        cont = 1
        while cont:
            shift = cell_nv[:,0] < regular_vertices.shape[0] 
            cell_nv[shift,:] = np.roll(cell_nv[shift,:],1,axis=1)
            cont = np.sum(shift) > 0
        shift = cell_nv[:,-1] < regular_vertices.shape[0] 
        cell_nv[shift,:] = np.roll(cell_nv[shift,:],-1,axis=1)
        
        
        cut_cells[nv-min_nv] = cell_nv
        cell_list_idx[nv-min_nv] = (cut_idx[0][idx], cut_idx[1][idx])
        irreg_edges[count:count+cell_nv.shape[0],:] = cell_nv[:,(0,-1)]
        count = count + cell_nv.shape[0]
    
    cell_list = [cells_whole] + cut_cells
    vertex_count = [4] + list(np.arange(min_nv,max_nv+1))
    cell_ij = [whole_idx] + cell_list_idx
    ncf = [0] + [1] * len(cut_cells)
    
    
    console.print("Computing the curved and corner boundaries...", style="bold blue")
    corners_flag = domain.is_corner(irreg_edges, vertices)
    num_corners = np.sum(corners_flag)
    regular_idx = np.where(np.logical_not(corners_flag))[0]
    corner_idx = np.where(corners_flag)[0]
    
    corner_cell_ij = [ () ] * len(cell_list)
    corner_cell_list = [np.array([])] * len(cell_list)
    corner_cut_nv = cut_nv[corner_idx]
    corner_vertex_count = [ v for v in vertex_count]
    # remove the corner cells
    for c in range(1,len(cell_list)):
        idx = np.where( cut_nv == vertex_count[c] )[0]
        
        keep = np.where(np.logical_not(corners_flag[idx]))[0]
        cidx = np.where(corners_flag[idx])[0]
    
        corner_cell_list[c] = cell_list[c][cidx, :]
        corner_cell_ij[c] = (cell_ij[c][0][cidx],cell_ij[c][1][cidx])
    
        
        cell_list[c] = cell_list[c][keep, :]
        cell_ij[c] = (cell_ij[c][0][keep],cell_ij[c][1][keep])

    cut_nv = cut_nv[regular_idx]
    

    if q > 1:
        # compute regular high order edges
        new_vertices, new_cells = domain.compute_curved(irreg_edges[regular_idx,:], vertices,q)
        shift =  vertices.shape[0]
        new_cells = new_cells + shift
        vertices = np.vstack( (vertices, new_vertices) )
        for c in range(1,len(cell_list)):
            idx = np.where( cut_nv == vertex_count[c] )[0]
            cell_list[c] = np.hstack( (cell_list[c], new_cells[idx,:]) ) 
            vertex_count[c] = vertex_count[c] + q-1




   
    if np.sum(corners_flag) > 0 :
        # deal with the corners
    
        
        new_vertices, new_cells = domain.compute_corner(irreg_edges[corner_idx,:], vertices,q)
        shift =  vertices.shape[0]
        new_cells = new_cells + shift
        vertices = np.vstack( (vertices, new_vertices) )
        
        count = 0
        for c in range(1,len(corner_cell_list)):
            if corner_cell_list[c].size == 0 :
                continue
            
            idx = corner_cell_list[c].shape[0]
            corner_cell_list[c] = np.hstack( (corner_cell_list[c], new_cells[count:count + idx,:]) ) 
            corner_vertex_count[c] = corner_cell_list[c].shape[1]
            count = count + idx
    
        # append corner bins
        cell_list = cell_list + corner_cell_list
        vertex_count = vertex_count + corner_vertex_count
        cell_ij = cell_ij + corner_cell_ij
        ncf = ncf + len(corner_cell_list) * [2] 


    # cell_list - bins of cells
    # vertex_count - number of vertices of cells in each bin
    # ncf - number of irregular faces associated to cells in each bin
    
    mesh_data = []
    for c, num_curved_faces in zip(cell_list, ncf):
        for idx in range(c.shape[0]):
            faces = []
            
            if num_curved_faces == 0:
                for vidx in range(c.shape[1]):
                    v1 = min([c[idx,vidx], c[idx,(vidx+1)%4]]) 
                    v2 = max([c[idx,vidx], c[idx,(vidx+1)%4]]) 
                    faces.append((v1,v2))
            else:
                num_extra = (q-1)*num_curved_faces+(num_curved_faces-1)
                final_vidx = c.shape[1]-num_extra

                for vidx in range(final_vidx-1):
                    v1=min([c[idx,vidx], c[idx,vidx+1]])
                    v2=max([c[idx,vidx], c[idx,vidx+1]])
                    faces.append( (v1, v2) )
                
                vertices_extra = list(c[idx,final_vidx-1:]) + [c[idx,0]]
                v1 = 0
                v2 =q+1
                for fidx in range(num_curved_faces):
                    faces.append( tuple(vertices_extra[v1:v2]) )
                    v1 = v2-1
                    v2 = v1 + q+1
            mesh_data.append({'vertices': c[idx,:], 'num_curved_faces': num_curved_faces, 'faces':faces})
    
    # order cells by number of faces
    face_number = [len(cell['faces']) for cell in mesh_data]
    elem_order = np.argsort(face_number)
    mesh_data = [mesh_data[idx] for idx in elem_order]
    for idx,elem in enumerate(mesh_data):
        elem['idx'] = idx
  
    face_hash = {}
    for elem in mesh_data:
        for f in elem['faces']:
            if f not in face_hash:
                face_hash[f] = {'lr':[elem['idx'], -1], 'vertices':f}
            else:
                face_hash[f]['lr'][1] = elem['idx']
    
    face_data = []
    for f in face_hash:
        face_data.append(face_hash[f])
    
    vertex_number = [len(face['vertices']) for face in face_data]
    face_order = np.argsort(vertex_number)
    face_data = [face_data[idx] for idx in face_order]
    for idx, f in enumerate(face_data):
        f['idx'] = idx

    for elem in mesh_data:
        elem['fidx'] = np.zeros((len(elem['faces']),)) 
        for idx, f in enumerate(elem['faces']):
            elem['fidx'][idx] = face_hash[f]['idx']

    
    
    
    
    
    
    # compute mesh stats
    def PolyArea(x,y):
        return 0.5*np.abs(np.sum(x * np.roll(y,1,axis=1),axis = 1)-np.sum(y * np.roll(x,1,axis=1), axis = 1))
    area = np.zeros((0,1))
    for cells,nv,ij in zip(cell_list, vertex_count, cell_ij):
        if cells.size == 0:
            continue
        
        xcoord = vertices[np.ravel(cells),0].reshape( (cells.shape[0], nv) )
        ycoord = vertices[np.ravel(cells),1].reshape( (cells.shape[0], nv) )
        
        in1 = xcoord >= XX[ij][:,None] 
        in2 = xcoord <= XX[ij[0]+1, ij[1]][:,None] 
        in3 = ycoord >= YY[ij][:,None] 
        in4 = ycoord <= YY[ij[0], ij[1]+1][:,None] 
    
        all_in = np.logical_and(in1,in2)
        all_in = np.logical_and(all_in,in3)
        all_in = np.logical_and(all_in, in4)
        any_out = np.logical_not(all_in)
        num_out = np.sum(any_out)
        
        if num_out > 0:
            print("There are some geomtrical vertices that lie outside a Cartesian cell")
            quit()
    
        ar = PolyArea(xcoord,ycoord)
        area = np.vstack( (area,ar[:,None]) )
        
    
    
    
    
    reg_vol = area[0]
    min_vol_frac = np.min(area) / reg_vol
    str_whole = str(cells_whole.shape[0])
    str_cut = str(num_cuts)
    str_tot = str(cells_whole.shape[0]+num_cuts)
    str_min_vol = '%.8e' % min_vol_frac[0]
    
    
    console.print("Embedded grid metadata", style="bold red")
    table = Table(show_header = False)
    table.add_column("Name", style="dim", width=15)
    table.add_column("Value",justify="right")
    table.add_row("whole cells", str_whole)
    table.add_row("cut cells", str_cut )
    table.add_row("total cells", str_tot )
    table.add_row("min. vol. frac.", str_min_vol )
    console.print(table)
    
    # reorder vertices counterclockwise
    if plot_flag:
        console.print("Plotting...", style="bold blue")
        plot_mesh(vertices,cell_list,vertex_count, dom, num_regular_vertices)
    
    return vertices, cell_list, domain, mesh_data, face_data
