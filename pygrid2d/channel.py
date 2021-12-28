import numpy as np
from .domain import Domain

class Channel(Domain):
    name = 'channel'
    left = 0
    right = 5.
    bottom = 0
    top = 5.
    num_blobs = 5


    def __init__(self):
        self.x1 = 6*5/30 + 0.5*5/30
        self.y1 = 0
        self.s1 = 1.

        self.x2 = 0
        self.y2 = 6*5/30 + 0.5*5/30
        self.s2 = 1.



    def in_domain(self, X,Y):
        Xp = X.ravel() 
        Yp = Y.ravel() 
        d1 =  -self.s1 * (Xp-self.x1) + (Yp-self.y1)
        d2 =  -self.s2 * (Xp-self.x2) + (Yp-self.y2)
        bval = np.logical_and( d1 > 0 , d2 < 0 ) 
        return bval.astype(int).reshape( X.shape )

    def bc(self,x,y):
        idx = -np.ones( (x.size) ).astype(int)
        return idx 
    def bc_id(self,bid):
        return bid[0]


    def vertex2bc(self, x, y):
        d1 =  -self.s1 * (x-self.x1) + (y-self.y1)
        d2 =  -self.s2 * (x-self.x2) + (y-self.y2)

        r1 = -1*(np.abs(d1)<1e-14) 
        r2 = -1*(np.abs(d2)<1e-14)
        left = -2*(np.abs(x-self.left)<1e-14)
        right = -2*(np.abs(x-self.right)<1e-14)
        top = -2*(np.abs(y-self.top)<1e-14)
        bot = -2*(np.abs(y-self.bottom)<1e-14)
        bc = r1+r2+top+bot+left+right
        
        decided = np.where(np.logical_xor(np.logical_xor(np.logical_xor(np.logical_xor(np.logical_xor(r1, r2), top), bot), left), right))
        out = np.zeros(x.shape)
        out[decided] = bc[decided]
        return out





    def curved_points(self,wgt,X1,Y1,X2,Y2):
        points = np.zeros( (X1.shape[0], wgt.size, 2) )
        points[:,:,0] = wgt[None,:]*X1[:,None] + (1.-wgt[None,:])*X2[:,None]
        points[:,:,1] = wgt[None,:]*Y1[:,None] + (1.-wgt[None,:])*Y2[:,None]
        return points 
