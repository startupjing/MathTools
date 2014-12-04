from __future__ import division
from pylab import *
from scipy.special.orthogonal import p_roots

def gauss1(f, n):
    """gaussian quadrature for weight=1 on (-1,1)"""
    [x,w] = p_roots(n+1)
    result = 0
    for i in range(n+1):
        result = result + w[i]*f(x[i])
    return result
        
def gauss(f, a, b, n):
    """gaussian quadrature for weight=1 on (a,b)"""
    [x,w] = p_roots(n+1)
    result = 0
    for i in range(n+1):
        temp = (b-a)*x[i]/2 + (b+a)/2
        result = result + w[i]*f(temp)
    return (b-a)*result/2
    
def f(x):
    return sqrt(1 - x**2)





