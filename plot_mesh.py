from matplotlib import collections  as mc
from matplotlib.collections import LineCollection
import numpy as np
import ipdb
import matplotlib.pyplot as plt



def plot_mesh(vertices, cells):
    ipdb.set_trace(context=21)
    segs = np.zeros((cells.shape[0], 5, 2))
    segs[:,0:4,0] = vertices[np.ravel(cells),0].reshape( (cells.shape[0], 4) )
    segs[:,0:4,1] = vertices[np.ravel(cells),1].reshape( (cells.shape[0], 4) )

    segs[:,4,:] = segs[:,0,:]

    line_segments = LineCollection(segs, linestyle='solid')
    
    ax = plt.axes()
    ax.set_xlim(vertices[:,0].min(), vertices[:,0].max())
    ax.set_ylim(vertices[:,1].min(), vertices[:,1].max())

    ax.add_collection(line_segments)
    plt.show()
