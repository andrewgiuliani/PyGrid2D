import numpy as np
from .domain import Domain

class Channel(Domain):
    name = 'channel'
    left = 0
    right = 5.
    bottom = 0
    top = 5.
    num_blobs = 5
    
#    def __init__(self):


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
        idx = -self.bc(X1,Y1) - 1

        radii = self.Ri[ 0, idx ]
        angle1 = np.arctan2( Y1-self.centroids[idx,1], X1-self.centroids[idx,0] )
        angle2 = np.arctan2( Y2-self.centroids[idx,1], X2-self.centroids[idx,0] )
        
       
        idx_jump = np.where(angle2-angle1 >   np.pi )
        angle1[idx_jump] = angle1[idx_jump] + 2.* np.pi
        idx_jump = np.where(angle2-angle1 <  -np.pi )
        angle2[idx_jump] = angle2[idx_jump] + 2.* np.pi


        angles = wgt[None,:]*angle1[:,None] + (1.-wgt[None,:])*angle2[:,None]        

        
        points = np.zeros( (angles.shape[0], angles.shape[1], 2) )
        points[:,:,0] = radii[:,None] * np.cos(angles)+self.centroids[idx,0][:,None]
        points[:,:,1] = radii[:,None] * np.sin(angles)+self.centroids[idx,1][:,None]
 
        return points 
        



