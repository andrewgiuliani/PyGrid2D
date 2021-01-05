import numpy as np
def fbody(X,Y,bid):

    body_ids = {
                0 :rh, 
                1 :ssv,
                2 :ringleb
                }
    active = body_ids[bid](X,Y)
    return active.astype(int)


def rh(X,Y):
    R1 = 0.75
    R2 = 1.25

    R = np.sqrt( (X-1.5)**2. + (Y-1.5)**2 )
    bval = np.logical_and( np.less_equal(R1 , R) , np.less_equal(R , R2) )
    return bval.astype(int)

def ssv(X,Y):
    R1 = 1.
    R2 = 1.384

    R = np.sqrt( X**2 + Y**2)
    bval = np.logical_and( np.less_equal(R1 , R) , np.less_equal(R , R2) )
    return bval.astype(int)

def ringleb(X,Y):
    return np.ones(X.shape)  
