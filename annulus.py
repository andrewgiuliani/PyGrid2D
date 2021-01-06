import numpy as np
from domain import Domain

class Annulus(Domain):
    
    name = 'annulus'
    left = 0
    right = 3.0001
    bottom = 0
    top = 3.0001

    
    def in_domain(self, X,Y):
        R1 = 0.75
        R2 = 1.25
    
        R = np.sqrt( (X-1.5)**2. + (Y-1.5)**2 )
        bval = np.logical_and( np.less_equal(R1 , R) , np.less_equal(R , R2) )
        return bval.astype(int)


