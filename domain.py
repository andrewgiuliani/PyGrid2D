import numpy as np
import ipdb
class Domain:

    def is_corner(self, irreg_edges, vertices):
        X1 = vertices[irreg_edges[:,0], 0]
        Y1 = vertices[irreg_edges[:,0], 1]
        X2 = vertices[irreg_edges[:,1], 0]
        Y2 = vertices[irreg_edges[:,1], 1]


        bc1 = self.bc(X1,Y1)
        bc2 = self.bc(X2,Y2)
        vals = np.logical_not(np.equal(bc1,bc2) )
        return vals
    
    def compute_curved(self, edges,vertices, q, bid):

        
        wgt = np.linspace(0.,1.,q+1)
        wgt = wgt[1:-1]
        X1 = vertices[edges[:,0],0]
        Y1 = vertices[edges[:,0],1]
        X2 = vertices[edges[:,1],0]
        Y2 = vertices[edges[:,1],1]
        extra_vertices=self.curved_points(wgt,X1,Y1,X2,Y2)
        
        out_vertices = extra_vertices.reshape( (-1, 2) )
        out_cells = np.arange(wgt.size * X1.size).reshape( (-1, wgt.size) )
        return out_vertices,out_cells
         
            
         
