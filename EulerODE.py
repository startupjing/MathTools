from __future__ import division
from pylab import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def euler(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        y[n+1] = y[n] + h*f(t[n], y[n])
    return y

def backwardEuler(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        def F(ynplus1):
            return -ynplus1 + y[n] + h*f(t[n+1], ynplus1)
        y[n+1] = fsolve(F, y[n])
    return y
    
def trapezoid(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        def F(ynplus1):
            return -ynplus1 + y[n] + 0.5*h*(f(t[n],y[n]) + f(t[n+1],ynplus1))
        y[n+1] = fsolve(F,y[n])
    return y

A = [[0,1],[-1,0]];

def f(t,y):
    return dot(A,y);
    
def main():
    y0 = [1,0];
    t0 = 0;
    tmax = 100;

    h1 = 0.1;
    N1 = 1000;
    y1 = euler(f,t0,y0,h1,N1);
    t1 = t0 + arange(N1+1)*h1;
    plt.plot(t1,y1[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('Euler Method with h=0.1');
    plt.show();
    
    h2 = 0.01;
    N2 = 10000;
    y2 = euler(f,t0,y0,h2,N2);
    t2 = t0 + arange(N2+1)*h2;
    plt.plot(t2,y2[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('Euler Method with h=0.01');
    plt.show();
    
    y3 = backwardEuler(f,t0,y0,h1,N1);
    plt.plot(t1,y3[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('BackwardEuler Method with h=0.1');
    plt.show();
    
    y4 = backwardEuler(f,t0,y0,h2,N2);
    plt.plot(t2,y4[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('BackwardEuler Method with h=0.01');
    plt.show();
    
    y5 = trapezoid(f,t0,y0,h1,N1);
    plt.plot(t1,y5[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('Trapezoid Method with h=0.1');
    plt.show();
    
    y6 = trapezoid(f,t0,y0,h2,N2);
    plt.plot(t2,y6[:,0]);
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('Trapezoid Method with h=0.01');
    plt.show();


if __name__ == "__main__": main()