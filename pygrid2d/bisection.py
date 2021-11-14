import numpy as np
def bisection(a,b,func):
    # you only need 55 iterations to get to machine precision
    max_step = 55
    
    init_diff = np.abs(a-b)

    for step in range(max_step):
        c = (a + b)/2.0
        fa = func(a)
        fc = func(c)
        fb = func(b)
    
        a = np.where( np.sign(fc) == np.sign(fa), c, a )
        b = np.where( np.sign(fc) == np.sign(fb), c, b )
    
    if init_diff.size > 0:
        max_res = np.max(init_diff)/2.**max_step
        if max_res > 1e-15:
            print("BISECTION NOT CONVERGED {}\n".format(max_res))
            quit()

    return c
