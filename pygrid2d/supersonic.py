import numpy as np
from .domain import Domain

class Supersonic(Domain):
    name = 'supersonic'
    left = 0
    right = 1.43
    bottom = 0
    top = 1.43
    R1 = 1.
    R2 = 1.384
 

    def in_domain(self, X,Y):
   
        R = np.sqrt( X**2 + Y**2)
        bval = np.logical_and( np.less_equal(self.R1 , R) , np.less_equal(R , self.R2) )
        return bval.astype(int)
    def bc_id(self,bid):
        return bid
    def bc(self,x,y):
        r = np.sqrt(x**2. +  y**2.)
        r1 = -1*(np.abs(r-self.R1) < 1e-14)
        r2 = -2*(np.abs(r-self.R2) < 1e-14)
        
        sanity_check = np.logical_not(np.logical_xor(r1,r2))
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("DUPLICATE BOUNDARY CONDITIONS\n")
            quit()
        
        return r1+r2    
    
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
        
