import numpy as np
import bisection as bs
from domain import Domain
import ipdb

class Ringleb(Domain):
    gamma = 1.4
    gm1 = gamma - 1.
    kmin = 0.7
    kmax = 1.2
    qmin = 0.5
    sx = 1.5
    sy = 0
    name = 'ringleb'    
    left = 0
    right = 2.75
    bottom = 0
    top = 2.75



    def compute_J(self,c):
        return 1. / c + 1. / (3. * c ** 3) + 1. / (5. * c ** 5) -.5 * np.log((1. + c) / (1. - c))
    def compute_rho(self,c):
        return c ** (2. / self.gm1) 
    def compute_q2(self,c):
        return 2. * (1. - c ** 2) / self.gm1
    
    def eqn(self,x, y, c):
        J = self.compute_J(c)
        rho = self.compute_rho(c)
        q2 = self.compute_q2(c)
    
        return (x - J / 2.) ** 2 + y ** 2 - 1. / (4. * rho ** 2 * q2 ** 2)
    
    def compute_c(self,x, y):
        clow =.5*np.ones(x.size)
        chigh = 0.9999*np.ones(x.size)
        
        c  = bs.bisection(clow, chigh, lambda cin : self.eqn(np.ravel(x),np.ravel(y),cin) )  
        return c.reshape(x.shape)
   
    def kq2xy(self, k, q):
        c = np.sqrt(2. + q**2 - self.gamma * q**2)/np.sqrt(2)
        J = self.compute_J(c)
        rho = self.compute_rho(c)
        q2 = q**2
        k2 = k**2
        
        x = (0.5/rho) * (1./q2 - 2./k**2.) + 0.5 * J + self.sx
        y = np.sqrt(1. - q2/k**2.)/(k * rho * q ) + self.sy
        return x,y
    
    def xy2kq(self,x,y):
        x = x - self.sx
        y = y - self.sy
    
        c = self.compute_c(x, y)
        J = self.compute_J(c)
        rho = self.compute_rho(c)
        q2 = self.compute_q2(c)
    
        q = np.sqrt(q2)
        k = np.sqrt(2. / (1. / q2 - 2. * rho * (x - J / 2.)))
        return k,q
    
    def in_domain(self,x,y):
        k,q = self.xy2kq(x,y)
        c1 = np.greater_equal(q,self.qmin)
        c2 = np.logical_and(np.less_equal(self.kmin, k), np.less_equal(k,self.kmax) )
        return np.logical_and(c1,c2).astype(int)

    def bc_id(self, bid):
        if bid == -1 or bid == -2:
            return -1
        else:
            return -2

    
    def bc(self,x,y):
        k,q = self.xy2kq(x,y)
        dk_min = -1*(np.abs(k-self.kmin) < 1e-13)
        dk_max = -2*(np.abs(k-self.kmax) < 1e-13)
        dq_min = -3*(np.abs(q-self.qmin) < 1e-13)
        
        sanity_check = np.logical_not(np.logical_xor(np.logical_xor(dk_min, dk_max), dq_min))
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("DUPLICATE BOUNDARY CONDITIONS\n")
            quit()
        
        return dk_min+dk_max+dq_min
    
    def target_ratio(self, wgt, k1, q1, k2, q2, k3, q3) :
        X1,Y1 =  self.kq2xy(k1,q1)
        X2,Y2 =  self.kq2xy(k2,q2)
        X3,Y3 =  self.kq2xy(k3,q3)
        
        l1 = np.sqrt( (X3-X2)**2. + (Y3-Y2)**2. )
        l3 = np.sqrt( (X3-X1)**2. + (Y3-Y1)**2. )
        
        return l1/l3 - wgt
        

    def curved_points(self,wgt,X1,Y1,X2,Y2):
        bc1 = self.bc(X1,Y1)
#        bc2 = self.bc(X2,Y2)
#        sanity_check = np.logical_not( np.equal(bc1,bc2) )
#        num_wrong = np.sum(sanity_check)
#        if(num_wrong > 0):
#            print("not the same\n")
#            quit()
        
        
        points = np.zeros( (X1.shape[0], wgt.size, 2) )
        for qq in range(wgt.size):
            # boundary -1
            idx1 = np.where( bc1 == -1 )[0]
            k1,q1 = self.xy2kq(X1[idx1], Y1[idx1])
            k3,q3 = self.xy2kq(X2[idx1], Y2[idx1])
            q2 = bs.bisection( q1, q3, lambda qin : self.target_ratio(wgt[qq], self.kmin, q1, 
                                                                               self.kmin, qin,
                                                                               self.kmin, q3 ) )
            x2,y2 = self.kq2xy(self.kmin,q2)
            points[idx1,qq,0] = x2
            points[idx1,qq,1] = y2

            # boundary -2
            idx2 = np.where( bc1 == -2)[0]
            k1,q1 = self.xy2kq(X1[idx2], Y1[idx2])
            k3,q3 = self.xy2kq(X2[idx2], Y2[idx2])
            q2 = bs.bisection( q1, q3, lambda qin : self.target_ratio(wgt[qq], self.kmax, q1, 
                                                                               self.kmax, qin,
                                                                               self.kmax, q3 ) )

            x2,y2 = self.kq2xy(self.kmax,q2)
            points[idx2,qq,0] = x2
            points[idx2,qq,1] = y2


            # boundary -3
            idx3 = np.where( bc1 == -3)[0]
            k1,q1 = self.xy2kq(X1[idx3], Y1[idx3])
            k3,q3 = self.xy2kq(X2[idx3], Y2[idx3])
            k2 = bs.bisection( k1, k3, lambda kin : self.target_ratio(wgt[qq], k1 , self.qmin, 
                                                                               kin, self.qmin,
                                                                               k3 , self.qmin ) )

            x2,y2 = self.kq2xy(k2,self.qmin)
            points[idx3,qq,0] = x2
            points[idx3,qq,1] = y2




        return points
 
    
    
    
    def get_corner_coord(self,bc1, bc2):
        k13,q13 = self.kq2xy(np.array([self.kmin]), np.array([self.qmin]))
        k23,q23 = self.kq2xy(np.array([self.kmax]), np.array([self.qmin]))
        
        corner13 = np.hstack( (k13,q13) ).reshape( (-1,2) )
        corner23 = np.hstack( (k23,q23) ).reshape( (-1,2) )

        corner_coords = np.zeros( (bc1.shape[0], 2) )
        corner_coords[:,0] = np.where( np.logical_and(bc1 == -1,bc2 == -3) , corner13[:,0], corner23[:,0] ) 
        corner_coords[:,1] = np.where( np.logical_and(bc1 == -1,bc2 == -3) , corner13[:,1], corner23[:,1] ) 

        return corner_coords





