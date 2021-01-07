import numpy as np
from domain import Domain
import ipdb

class Blobs(Domain):
    name = 'blobs'
    left = 0
    right = 1.
    bottom = 0
    top = 1.
    num_blobs = 10 
    num_modes = 4
    
    def __init__(self):
        #XX,YY = np.meshgrid( np.linspace(self.left,self.right,self.num_blobs) , np.linspace(self.bottom,self.top, self.num_blobs))
        #self.centroids = np.hstack( (XX.ravel()[:,None], YY.ravel()[:,None] ) ) 

        

        self.centroids = 0.1 + 0.9 * np.random.rand( self.num_blobs,2  )
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
        
        bval = np.sum(R <= R_blob, axis = 1) > 0

#        ipdb.set_trace(context=21)
        return bval.astype(int).reshape( X.shape )

    def bc(self,x,y):
        return -np.ones(x.shape) 
    
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        bc1 = self.bc(X1,Y1)
        bc2 = self.bc(X2,Y2)
        sanity_check = np.logical_not( np.equal(bc1,bc2) )
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("not the same\n")
            quit()
        
        radii = np.where(np.equal(bc1, -1) , self.R1, self.R2)
        angle1 = np.arctan2( Y1, X1 )
        angle2 = np.arctan2( Y2, X2 )
        angles = wgt[None,:]*angle1[:,None] + (1.-wgt[None,:])*angle2[:,None]        
        
        points = np.zeros( (angles.shape[0], angles.shape[1], 2) )
        points[:,:,0] = radii[:,None] * np.cos(angles)
        points[:,:,1] = radii[:,None] * np.sin(angles)
        return points
        



