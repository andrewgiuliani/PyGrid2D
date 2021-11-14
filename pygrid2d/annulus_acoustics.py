import numpy as np
from .domain import Domain

class Annulus_acoustics(Domain):
    
    name = 'annulus_acoustics'
    left = 0.
    right = 21.
    bottom = 0.
    top = 21.
    R1 = 0.5
    R2 = 10.
    sx = 10.5
    sy = 10.5
 
    
    def in_domain(self, X,Y):
        R = np.sqrt( (X-self.sx)**2. + (Y-self.sy)**2 )
        bval = np.logical_and( np.less_equal(self.R1 , R) , np.less_equal(R , self.R2) )
        return bval.astype(int)

    def bc(self,x,y):
        r = np.sqrt( (x-self.sx)**2. +  (y-self.sy)**2.)
        r1 = -1*(np.abs(r-self.R1) < 1e-14)
        r2 = -2*(np.abs(r-self.R2) < 1e-14)
        
        sanity_check = np.logical_not(np.logical_xor(r1,r2))
        num_wrong = np.sum(sanity_check)
        
        if(num_wrong > 0):
            print("DUPLICATE BOUNDARY CONDITIONS\n")
            quit()
        
        return r1+r2    

    def bc_id(self, bid):
        return bid
    
    def curved_points(self,wgt,X1,Y1,X2,Y2):
        bc1 = self.bc(X1,Y1)
        bc2 = self.bc(X2,Y2)
        sanity_check = np.logical_not( np.equal(bc1,bc2) )
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("not the same\n")
            quit()
        
        radii = np.where(np.equal(bc1, -1) , self.R1, self.R2)
        angle1 = np.arctan2( Y1-self.sy, X1-self.sx )
        angle2 = np.arctan2( Y2-self.sy, X2-self.sx )
        
       
        idx_jump = np.where(angle2-angle1 >   np.pi )
        angle1[idx_jump] = angle1[idx_jump] + 2.* np.pi
        idx_jump = np.where(angle2-angle1 <  -np.pi )
        angle2[idx_jump] = angle2[idx_jump] + 2.* np.pi


        angles = wgt[None,:]*angle1[:,None] + (1.-wgt[None,:])*angle2[:,None]        

        
        points = np.zeros( (angles.shape[0], angles.shape[1], 2) )
        points[:,:,0] = radii[:,None] * np.cos(angles)+self.sx
        points[:,:,1] = radii[:,None] * np.sin(angles)+self.sy
        return points
 
