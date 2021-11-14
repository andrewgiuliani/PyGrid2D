import numpy as np
from .domain import Domain

class Channel(Domain):
    name = 'channel'
    left = 0
    right = 5.
    bottom = 0
    top = 5.
    num_blobs = 5
    
    def in_domain(self, X,Y):
        Xp = X.ravel() 
        Yp = Y.ravel() 
        x1 = 6*5/30 + 0.5*5/30
        y1 = 0
        s1 = 1.

        x2 = 0
        y2 = 6*5/30 + 0.5*5/30
        s2 = 1.

        d1 =  -s1 * (Xp-x1) + (Yp-y1)
        d2 =  -s2 * (Xp-x2) + (Yp-y2)
        bval = np.logical_and( d1 > 0 , d2 < 0 ) 
        return bval.astype(int).reshape( X.shape )

    def bc(self,x,y):
        idx = -np.ones( (x.size) ).astype(int)
        return idx 
    def bc_id(self,bid):
        return bid[0]
    
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        points = np.zeros( (X1.shape[0], wgt.size, 2) )
        points[:,:,0] = wgt[None,:]*X1[:,None] + (1.-wgt[None,:])*X2[:,None]
        points[:,:,1] = wgt[None,:]*Y1[:,None] + (1.-wgt[None,:])*Y2[:,None]
        return points 
