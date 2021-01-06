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
   
    def boundary_curve(self, k, q):
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
    
    
    def bc(self,x,y):
        k,q = self.xy2kq(x,y)
        dk_min = 1*(np.abs(k-kmin) < 1e-14)
        dk_max = 2*(np.abs(k-kmax) < 1e-14)
        dq_min = 3*(np.abs(q-qmin) < 1e-14)
        
        sanity_check = np.logical_not(np.logical_xor(np.logical_xor(dk_min, dk_max), dq_min))
        num_wrong = np.sum(sanity_check)
        if(num_wrong > 0):
            print("DUPLICATE BOUNDARY CONDITIONS\n")
            quit()
        
        return dk_min+dk_max+dq_min
    
    
    
    
       
