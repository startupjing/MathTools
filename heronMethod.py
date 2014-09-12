from __future__ import division
from pylab import *

def heronStep(y, x):
    """Perform one step of Heron's method for sqrt(y) with guess x."""
    return 0.5*(x + y/x)

def heronArray(y, x0, n):
    """Perform n steps of Heron's method for sqrt(y) with initial guess x0.
    Return the array of successive approximations [x0, x1, ..., xn]."""
    x = zeros(n+1)
    x[0] = x0
    for k in range(n):
        x[k+1] = heronStep(y, x[k])
    return x

def cosArray(x0, n):
    x = zeros(n+1)
    x[0] = x0
    for k in range(n):
        x[k+1] = cos(x[k])
    return x

def plotLogError(y, x0, n):
    """Plot the log absolute error of Heron's method."""
    x = heronArray(y, x0, n) # Heron's method approximations for sqrt(y)
    e = x - sqrt(y) # error
    semilogy(abs(e)) # plot abs(e) using a log scale for the y-axis
