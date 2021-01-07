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
    
    def compute_curved(self, edges,vertices, q):

        
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
         
            
    def compute_corner(self, edges, vertices, q):
        X1 = vertices[edges[:,0], 0]
        Y1 = vertices[edges[:,0], 1]
        X2 = vertices[edges[:,1], 0]
        Y2 = vertices[edges[:,1], 1]
        bc1 = self.bc(X1,Y1)
        bc2 = self.bc(X2,Y2)
        sanity_check = np.equal(bc1,bc2)
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("not the same\n")
            quit()
        num_corners = edges.shape[0]
        new_corner_coord = self.get_corner_coord(bc1,bc2)
        
        wgt = np.linspace(0.,1.,q+1)
        wgt1 = wgt[:-1]
        wgt2 = wgt[1:-1]
        extra_vertices1=self.curved_points(wgt1, X1,Y1,new_corner_coord[:,0],new_corner_coord[:,1])
        extra_vertices2=self.curved_points(wgt2, X2,Y2,new_corner_coord[:,0],new_corner_coord[:,1])
        extra_vertices1[:,:,0] = np.fliplr(extra_vertices1[:,:,0])
        extra_vertices1[:,:,1] = np.fliplr(extra_vertices1[:,:,1])
#        extra_vertices2[:,:,0] = np.fliplr(extra_vertices2[:,:,0])
#        extra_vertices2[:,:,1] = np.fliplr(extra_vertices2[:,:,1])
       
        extra_vertices1 = extra_vertices1.reshape( (-1,2) )
        extra_vertices2 = extra_vertices2.reshape( (-1,2) )
        extra_vertices = np.vstack( (extra_vertices1, extra_vertices2) )

        extra_cells1 = np.arange( 0, X1.shape[0] * q ).reshape( (-1, q ) )
        extra_cells2 = np.zeros( (extra_cells1.shape[0],0) ).astype(int)
        if q > 1:
            extra_cells2 = np.arange( 0, X1.shape[0] * (q-1) ).reshape( (-1, q-1 ) ) + np.max(extra_cells1)+1
        extra_cells = np.fliplr(np.hstack( ( extra_cells1, extra_cells2) ))
        
        return extra_vertices, extra_cells
