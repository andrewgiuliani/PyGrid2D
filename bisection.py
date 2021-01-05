import numpy as np
def bisection(a,b,func):
    max_step = 75
    

    for step in range(max_step):
        c = (a + b)/2.0
        fa = func(a)-0.5
        fc = func(c)-0.5
        fb = func(b)-0.5
    
        a = np.where( np.sign(fc) == np.sign(fa), c, a )
        b = np.where( np.sign(fc) == np.sign(fb), c, b )

    return c
