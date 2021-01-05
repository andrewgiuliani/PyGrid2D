import numpy as np
def in_domain(X,Y):
    R1 = 1.
    R2 = 1.384

    R = np.sqrt( X**2 + Y**2)
    bval = np.logical_and( np.less_equal(R1 , R) , np.less_equal(R , R2) )
    return bval.astype(int)


