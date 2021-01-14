import numpy as np
from domain import Domain
import bisection as bs
import ipdb

class Blobs(Domain):
    name = 'blobs'
    left = 0
    right = 1.
    bottom = 0
    top = 1.
    num_blobsx = 3
    num_blobsy = 3
    num_modes = 3
    
    def __init__(self):
        XX,YY = np.meshgrid( np.linspace(self.left,self.right,self.num_blobsx) , np.linspace(self.bottom,self.top, self.num_blobsy))
        XX = 0.6 * XX + 0.15
        YY = 0.6 * YY + 0.15
        self.centroids = np.hstack( (XX.ravel()[:,None], YY.ravel()[:,None] ) ) 
        self.num_blobs = self.centroids.shape[0]
        

 #       ipdb.set_trace(context=21)
#        self.centroids = 0.1 + 0.9 * np.random.rand( self.num_blobs,2  )
        self.modes_sin = np.random.rand( self.num_blobs, self.num_modes )
        self.modes_cos = np.random.rand( self.num_blobs, self.num_modes )
        
        scale = np.array(  range(1,self.num_modes+1) ) * 10
        scale[0] = 1
        self.modes_sin[:,0] = 0.1
        self.modes_cos[:,0] = 0.1
        self.modes_sin = self.modes_sin / scale
        self.modes_cos = self.modes_cos / scale


    def in_domain(self, X,Y):
        Xp = X.ravel() 
        Yp = Y.ravel() 
        R = np.sqrt( (Xp[:,None] - self.centroids[:,0])**2 + (Yp[:,None] - self.centroids[:,1])**2)
        theta = np.arctan2(Yp[:,None] - self.centroids[:,1], Xp[:,None] - self.centroids[:,0])
        
        ntheta = np.array( range(self.num_modes) )[None,None,:] * theta[:,:,None]
        R_blob =  np.sum(self.modes_sin[None,:,:] * np.sin(ntheta) + self.modes_cos[None,:,:] * np.cos(ntheta), axis = 2)
        
        bval = np.logical_not(np.sum(R <= R_blob, axis = 1) > 0)
        return bval.astype(int).reshape( X.shape )

    def bc(self,x,y):
        Xp = x.ravel() 
        Yp = y.ravel() 
        R = np.sqrt( (Xp[:,None] - self.centroids[:,0])**2 + (Yp[:,None] - self.centroids[:,1])**2)
        theta = np.arctan2(Yp[:,None] - self.centroids[:,1], Xp[:,None] - self.centroids[:,0])
        
 
        ntheta = np.array( range(self.num_modes) )[None,None,:] * theta[:,:,None]
        R_blob =  np.sum(self.modes_sin[None,:,:] * np.sin(ntheta) + self.modes_cos[None,:,:] * np.cos(ntheta), axis = 2)
 
        rdiff = np.sqrt( (Xp[:,None]-self.centroids[:,0])**2. +  (Yp[:,None]-self.centroids[:,1])**2.)-R_blob
        bools = rdiff < 1e-14
        idx = -np.where(rdiff < 1e-14)[1]-1
        return idx 
    def bc_id(self,bid):
        return -np.ones(bid.shape)
  

    def target_ratio(self, wgt, idx,R1, angle1, angle2, R3, angle3) :
        X1,Y1 = R1 * np.cos(angle1), R1* np.sin(angle1) 
        
        ntheta = np.array( range(self.num_modes) )[None,:] * angle2[:,None]
        R2 =  np.sum(self.modes_sin[idx,:] * np.sin(ntheta) + self.modes_cos[idx,:] * np.cos(ntheta), axis = 1)
        X2,Y2 = R2 * np.cos(angle2), R2 * np.sin(angle2)  
        
        X3,Y3 = R3 * np.cos(angle3), R3 * np.sin(angle3) 
 
        
        l1 = np.sqrt( (X3-X2)**2. + (Y3-Y2)**2. )
        l3 = np.sqrt( (X3-X1)**2. + (Y3-Y1)**2. )
        
        return l1/l3 - wgt
 
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        idx = -self.bc(X1,Y1) - 1

        R1 = np.sqrt( (X1 - self.centroids[idx,0])**2 + (Y1 - self.centroids[idx,1])**2)
        R2 = np.sqrt( (X2 - self.centroids[idx,0])**2 + (Y2 - self.centroids[idx,1])**2)
        angle1 = np.arctan2(Y1 - self.centroids[idx,1], X1 - self.centroids[idx,0])
        angle2 = np.arctan2(Y2 - self.centroids[idx,1], X2 - self.centroids[idx,0])
          
       
        idx_jump = np.where(angle2-angle1 >   np.pi )
        angle1[idx_jump] = angle1[idx_jump] + 2.* np.pi
        idx_jump = np.where(angle2-angle1 <  -np.pi )
        angle2[idx_jump] = angle2[idx_jump] + 2.* np.pi

        
        points = np.zeros( (X1.shape[0], wgt.size, 2) )
        for qq in range(wgt.size):
            angle3 = bs.bisection( angle1, angle2, lambda qin : self.target_ratio(wgt[qq], idx,R1, angle1, qin, R2, angle2) )
#            ipdb.set_trace(context=21)
            
            ntheta = np.array( range(self.num_modes) )[None,:] * angle3[:,None]
            R3 =  np.sum(self.modes_sin[idx,:] * np.sin(ntheta) + self.modes_cos[idx,:] * np.cos(ntheta), axis = 1)
            X3,Y3 = R3 * np.cos(angle3)+self.centroids[idx,0], R3 * np.sin(angle3) + self.centroids[idx,1] 
             
            points[:,qq,0] = X3
            points[:,qq,1] = Y3



        return points 
 



#    def curved_points(self,wgt,X1,Y1,X2,Y2):
#        ipdb.set_trace(context=21)
#        idx = -self.bc(X1,Y1) - 1
#
#        radii = self.Ri[ 0, idx ]
#        angle1 = np.arctan2( Y1-self.centroids[idx,1], X1-self.centroids[idx,0] )
#        angle3 = np.arctan2( Y2-self.centroids[idx,1], X2-self.centroids[idx,0] )
#        
#       
#        idx_jump = np.where(angle3-angle1 >   np.pi )
#        angle1[idx_jump] = angle1[idx_jump] + 2.* np.pi
#        idx_jump = np.where(angle3-angle1 <  -np.pi )
#        angle3[idx_jump] = angle3[idx_jump] + 2.* np.pi
#
#
 
#        return points 
#        
#
#
#
