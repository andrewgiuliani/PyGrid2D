import numpy as np
from .domain import Domain

class DoubleMach(Domain):
    name = 'doublemach'
    left = 0
    right = 2.5
    bottom = 0
    top = 1.75
    corner_x = 1./6.
    corner_y = 1e-10

    def in_domain(self, X,Y):
        dist_from_wedge = (X - self.corner_x) * np.sin(30*np.pi/180) + (Y - self.corner_y) * -np.cos(30 * np.pi/180)
        bval1 = np.less_equal(dist_from_wedge, 0)
        bval2 = np.less_equal(self.corner_y, Y)
        bval = np.logical_and(bval1, bval2)
        return bval.astype(int)
    
    def bc_id(self,bid):
        return bid
    
    def bc(self,x,y):
        dist_from_wedge = (x - self.corner_x) * np.sin(30*np.pi/180) + (y - self.corner_y) * -np.cos(30 * np.pi/180)
        bc_out = -1*(np.abs(y-self.corner_y) < 1e-13) - 2 * (np.abs(dist_from_wedge) < 1e-13)
        return bc_out.astype(int)
    
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        points = np.zeros( (X1.shape[0], wgt.size, 2) )
        points[:,:,0] = wgt[None,:]*X1[:,None] + (1.-wgt[None,:])*X2[:,None]
        points[:,:,1] = wgt[None,:]*Y1[:,None] + (1.-wgt[None,:])*Y2[:,None]
        return points
    
    def get_corner_coord(self,bc1, bc2):
        return np.array([self.corner_x, self.corner_y]).reshape((-1,2))

