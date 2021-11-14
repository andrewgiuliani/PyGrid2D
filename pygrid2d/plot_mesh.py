from matplotlib import collections  as mc
from matplotlib.collections import LineCollection
import numpy as np
import matplotlib.pyplot as plt



def plot_mesh(vertices, cell_list, vertex_counts, dom, regular_vertices):
    
    ax = plt.axes()
    ax.set_xlim(dom[0], dom[2])
    ax.set_ylim(dom[1], dom[3])
    
#    cell_list = [cell_list[1]]
#    vertex_counts = [vertex_counts[1]]
    for cells,nv in zip(cell_list,vertex_counts):
        if cells.size == 0:
            continue

        #ipdb.set_trace(context=21)
        segs = np.zeros((cells.shape[0], nv+1, 2))
        segs[:,0:nv,0] = vertices[np.ravel(cells),0].reshape( (cells.shape[0], nv) )
        segs[:,0:nv,1] = vertices[np.ravel(cells),1].reshape( (cells.shape[0], nv) )

        segs[:,nv,:] = segs[:,0,:]
        line_segments = LineCollection(segs, linestyle='solid', linewidth = 0.5)
        
        
        ax.add_collection(line_segments)
    ax.set_aspect('equal', 'box')
    
    plt.scatter(vertices[regular_vertices:,0], vertices[regular_vertices:,1], color = 'red', s = 1) 
    plt.scatter(vertices[0:regular_vertices,0], vertices[0:regular_vertices,1], color = 'black', s = 1) 
    plt.show()
