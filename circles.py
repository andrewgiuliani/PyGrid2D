import numpy as np
from domain import Domain
import ipdb

class Circles(Domain):
    name = 'circles'
    left = 0
    right = 20.
    bottom = 0
    top = 20.
    sx = 10
    sy = 10
    num_blobs = 5
    
    def __init__(self):
        theta = np.linspace( 0, 2 * np.pi, self.num_blobs+1)[:,None]
        theta = theta[:-1]
        Ro = 6
        self.centroids = np.hstack( (Ro * np.cos(theta)+self.sx, Ro*np.sin(theta)+self.sy ) ) 
        self.num_blobs = self.centroids.shape[0]
        self.Ri = 2*np.ones ((1,self.num_blobs) )


    def in_domain(self, X,Y):
        Xp = X.ravel() 
        Yp = Y.ravel() 
#        ipdb.set_trace(context=21)
        R = np.sqrt( (Xp[:,None] - self.centroids[None,:,0])**2 + (Yp[:,None] - self.centroids[None,:,1])**2)
        theta = np.arctan2(Yp[:,None] - self.centroids[:,1], Xp[:,None] - self.centroids[:,0])
        
        bval = np.logical_not(np.sum(R <= self.Ri, axis = 1) > 0)
        return bval.astype(int).reshape( X.shape )

    def bc(self,x,y):
        x = np.array(x).reshape( (-1,1) )
        y = np.array(y).reshape( (-1,1) )

        rdiff = np.sqrt( (x-self.centroids[:,0])**2. +  (y-self.centroids[:,1])**2.)-self.Ri
        bools = rdiff < 1e-14
        idx = -np.where(rdiff < 1e-14)[1]-1
        return idx 
    def bc_id(self,bid):
        return bid
    
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
 
        #ipdb.set_trace(context=21)
        return points 
        



