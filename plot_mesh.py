from matplotlib import collections  as mc
from matplotlib.collections import LineCollection
import numpy as np
import ipdb
import matplotlib.pyplot as plt



def plot_mesh(vertices, cell_list, vertex_counts):
    
    ax = plt.axes()
    ax.set_xlim(vertices[:,0].min(), vertices[:,0].max())
    ax.set_ylim(vertices[:,1].min(), vertices[:,1].max())
    
#    cell_list = [cell_list[1]]
#    vertex_counts = [vertex_counts[1]]
    for cells,nv in zip(cell_list,vertex_counts):
        #ipdb.set_trace(context=21)
        
        segs = np.zeros((cells.shape[0], nv+1, 2))
        segs[:,0:nv,0] = vertices[np.ravel(cells),0].reshape( (cells.shape[0], nv) )
        segs[:,0:nv,1] = vertices[np.ravel(cells),1].reshape( (cells.shape[0], nv) )

        segs[:,nv,:] = segs[:,0,:]

        line_segments = LineCollection(segs, linestyle='solid')
        
        
        ax.add_collection(line_segments)
    plt.show()
