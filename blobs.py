import numpy as np
from domain import Domain
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
        XX = 0.6 * XX + 0.2
        YY = 0.6 * YY + 0.2
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
        return -np.ones(x.shape) 
    
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        print("No curved points implemented yet")
        quit()




