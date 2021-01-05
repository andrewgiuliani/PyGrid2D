import numpy as np
import bisection as bs
import ipdb

gamma = 1.4
gm1 = gamma - 1.
kmin = 0.7
kmax = 1.2
qmin = 0.5
sx = 1.5
sy = 0

def compute_J(c):
    return 1. / c + 1. / (3. * c ** 3) + 1. / (5. * c ** 5) -.5 * np.log((1. + c) / (1. - c))
def compute_rho(c):
    return c ** (2. / gm1) 
def compute_q2(c):
    return 2. * (1. - c ** 2) / gm1

def eqn(x, y, c):
    J = compute_J(c)
    rho = compute_rho(c)
    q2 = compute_q2(c)

    return (x - J / 2.) ** 2 + y ** 2 - 1. / (4. * rho ** 2 * q2 ** 2)

def compute_c(x, y):
    clow =.5*np.ones(x.size)
    chigh = 0.9999*np.ones(x.size)
    
    c  = bs.bisection(clow, chigh, lambda cin : eqn(np.ravel(x),np.ravel(y),cin) )  
    return c.reshape(x.shape)


def xy2kq(x,y):
    x = x - sx
    y = y - sy

    c = compute_c(x, y)
    J = compute_J(c)
    rho = compute_rho(c)
    q2 = compute_q2(c)

    q = np.sqrt(q2)
    k = np.sqrt(2. / (1. / q2 - 2. * rho * (x - J / 2.)))
    return k,q

def in_domain(x,y):
    k,q = xy2kq(x,y)
    c1 = np.greater_equal(q,qmin)
    c2 = np.logical_and(np.less_equal(kmin, k), np.less_equal(k,kmax) )
    return np.logical_and(c1,c2)


def bc(k,q):
    dk_min = 1*(np.abs(k-kmin) < 1e-13)
    dk_max = 2*(np.abs(k-kmax) < 1e-13)
    dq_min = 3*(np.abs(q-qmin) < 1e-13)
    
    sanity_check = np.logical_not(np.logical_xor(np.logical_xor(dk_min, dk_max), dq_min))
    num_wrong = np.sum(sanity_check)
    if(num_wrong > 0):
        print("DUPLICATE BOUNDARY CONDITIONS\n")
        quit()
    
    return dk_min+dk_max+dq_min




def is_corner(X1,Y1, X2,Y2):
    k1,q1 = xy2kq(X1,Y1)
    k2,q2 = xy2kq(X2,Y2)

    bc1 = bc(k1,q1)
    bc2 = bc(k2,q2)
    
    vals = np.logical_not(np.equal(bc1,bc2) )
    return vals
    
