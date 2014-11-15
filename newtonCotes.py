from __future__ import division
from pylab import *

def mid(f, a, b):
    """Midpoint Rule for function f on [a,b]"""
    return 1.0*(b-a)*f((a+b)/2.0)
    
def trap(f, a, b):
    """Trapezoid Rule for function f on [a,b]"""
    return 0.5*(b-a)*(f(a) + f(b))
    
def simp(f, a, b):
    """Simpson's Rule for function f on [a,b]"""
    return (b-a)*(f(a) + 4*f((a+b)/2.0) + f(b))/6.0

def midc(f, a, b, m):
    """Composite midpoint rule for function on [a,b]"""
    h = 1.0*(b-a)/m
    sum = 0
    for i in range(1,m+1):
        x1 = a + (i-1)*h
        x2 = a + i*h
        sum += f((x1+x2)/2.0)
    return h*sum
    
def trapc(f, a, b, m):
    """Composite trapezoid rule for function f on [a,b]"""
    h = 1.0*(b-a)/m
    sum = 0
    for i in range(1,m+1):
        x = a + i*h
        if i==0 or i==m:
            sum += 0.5*f(x)
        else:
            sum += f(x)
    return h*sum
    
def simpc(f, a, b, m):
    """Composite simpson's rule for function f on [a,b]"""
    h = 0.5*(b-a)/m
    sum = 0
    for i in range(1,m+1):
        x1 = a + (2*i -2)*h
        x2 = a + (2*i -1)*h
        x3 = a + 2*i*h
        sum += x1+x2+x3
    return sum*h/3
        
def f(x):
    return sqrt(1-x**2)
    
    

