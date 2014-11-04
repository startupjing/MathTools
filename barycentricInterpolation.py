from __future__ import division
from pylab import *

def weights(xk):
    """Function takes the array of interpolation pointx
       xk = (x0, x1, ..., xn) as input and returns the array
       of barycentric weights wk = (w0, w1, ..., wn) as output"""
    wk = ones(xk.size)
    for k in range(wk.size):
        for i in range(wk.size):
            if i != k:
                wk[k] = wk[k] * (1.0/(xk[k] - xk[i]));
    return wk
    
def interpolate(x, xk, yk, wk):
    """Function takes x as the evaluation point, xk as the array
       of interpolation points (x0,x1,...,xn), yk as the array of
       interpolation values (y0,y1,...,yn), and wk as the array of
       barycentric weights (w0,w1,...,wn). The function returns the
       value of pn(x)"""
    if x in xk:
        return null
    else:
        numer = 0;
        denom = 0;
        temp = 0;
        for k in range(xk.size):
            temp = 1.0* wk[k] / (x - xk[k])
            numer = numer + temp*yk[k]
            denom = denom + temp
        return numer/denom
        
def f(x):
    """Define function f(x)= 1/(1+x^2)"""
    return 1.0/(1+x**2)
            
    
    

