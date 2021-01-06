import numpy as np
def bisection(a,b,func):
    # you only need 55 iterations to get to machine precision
    max_step = 55
    

    for step in range(max_step):
        c = (a + b)/2.0
        fa = func(a)
        fc = func(c)
        fb = func(b)
    
        a = np.where( np.sign(fc) == np.sign(fa), c, a )
        b = np.where( np.sign(fc) == np.sign(fb), c, b )

    return c
