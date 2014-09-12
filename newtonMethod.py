from __future__ import division
from pylab import *

def newtonStep(f, df, x):
    """Perform one step of Newton's method."""
    return x - f(x)/df(x)

def newtonArray(f, df, x0, n):
    """Perform n steps of Newton's method with initial guess x0.
    Return the array of successive approximations [x0, x1, ..., xn]."""
    x = zeros(n+1)
    x[0] = x0
    for k in range(n):
        x[k+1] = newtonStep(f, df, x[k])
    return x

def recipStep(y, x):
    """Perform one step of reciprocal algorithm."""
    return 2*x - y*(x**2);


def recipArray(y, x0, n):
    """Perform n steps of reciprocal algorithm with initial guess x0.
       Return the array of successive approximations [x0, x1, ..., xn]"""
    x = zeros(n+1)
    x[0] = x0
    for k in range(n):
        x[k+1] = recipStep(y, x[k])
    return x

def f(x):
    return x**3 - 2
    
def df(x):
    return 3*(x**2)