from __future__ import division
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def euler(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        y[n+1] = y[n] + h*f(t[n], y[n])
    return y

def etrap(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        xi1 = y[n]
        f1 = f(t[n], xi1)
        xi2 = y[n] + h*f1
        f2 = f(t[n+1], xi2)
        y[n+1] = y[n] + 0.5*h*(f1 + f2)
    return y
    
def emid(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0]= y0
    for n in range(N):
        xi1 = y[n]
        f1 = 0
        xi2 = y[n] + 0.5*h*f(t[n+1], xi1)
        f2 = f(t[n]+0.5*h, xi2)
        y[n+1] = y[n] + h*f2
    return y
    
def rk3(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        xi1 = y[n]
        f1 = f(t[n], xi1)
        xi2 = y[n] + 0.5*h*f(t[n], xi1)
        f2 = f(t[n]+0.5*h, xi2)
        xi3 = y[n] - h*f(t[n], xi1) + 2*h*f(t[n]+0.5*h, xi2)
        f3 = f(t[n+1], xi3)
        y[n+1] = y[n] + h*(f1/6 + 2*f2/3 + f3/6)
    return y
    

def rk4(f, t0, y0, h, N):
    t = t0 + arange(N+1)*h
    y = zeros((N+1, size(y0)))
    y[0] = y0
    for n in range(N):
        xi1 = y[n]
        f1 = f(t[n], xi1)
        xi2 = y[n] + 0.5*h*f(t[n], xi1)
        f2 = f(t[n]+0.5*h, xi2)
        xi3 = y[n] + 0.5*h*f(t[n]+0.5*h, xi2)
        f3 = f(t[n]+0.5*h, xi3)
        xi4 = y[n] + h*f(t[n]+0.5*h, xi3)
        f4 = f(t[n+1], xi4)
        y[n+1] = y[n] + h*(f1/6 + f2/3 + f3/3 + f4/6)
    return y
    
def f(t,y):
    return y

def errorPlot(method):
    def f(t,y):
        return y
    err = zeros(10)
    H = 0.5**arange(10)
    for i, h in enumerate(H):
        N = int(1/h)
        y = method(f, 0, 1, h, N)
        err[i] = abs(y[N] - e)
    loglog(H,err)
    xlabel("h")
    ylabel("error")
    
def plots(method1, method2, method3, method4):
    def f(t,y):
        return y
    err1 = zeros(10)
    err2 = zeros(10)
    err3 = zeros(10)
    err4 = zeros(10)
    H = 0.5**arange(10)
    for i, h in enumerate(H):
        N = int(1/h)
        y1 = method1(f, 0, 1, h, N)
        err1[i] = abs(y1[N] - e)
        y2 = method2(f, 0, 1, h, N)
        err2[i] = abs(y2[N] - e)
        y3 = method3(f, 0, 1, h, N)
        err3[i] = abs(y3[N] - e)
        y4 = method4(f, 0, 1, h, N)
        err4[i] = abs(y4[N] - e)
    plt.loglog(H,err1,'g--',label='euler')
    plt.loglog(H,err2,'r-o',label='emid')
    plt.loglog(H,err3,'+',label='rk3')
    plt.loglog(H,err4,'x',label='rk4')
    xlabel("h")
    ylabel("error")
    plt.legend()
    plt.show()


def fLorenz(t, y):
    val = zeros(size(y))
    val[0] = 10*(y[1] - y[0])
    val[1] = y[0]*(28 - y[2]) - y[1]
    val[2] = y[0]*y[1] - 8*y[2]/3
    return val
    
    
def lorenzPlot():
    y = rk4(fLorenz, 0, array([0,2,20]), .01, 10000)
    fig = figure()
    ax = Axes3D(fig)
    ax.plot(*y.T)