from __future__ import division
from pylab import *
from scipy.linalg import solveh_banded

b = [1]*10000

diag = [2]*10000
sub_diag = [-1]*9999
ab = array([[0]+sub_diag, diag])

def tridiag(sub_diag1, diag, sub_diag2, k1=-1, k2=0, k3=1):
    return np.diag(sub_diag1, k1) + np.diag(diag, k2) + np.diag(sub_diag2, k3)

a = tridiag(sub_diag, diag, sub_diag)


#x_h = solveh_banded(ab,b)
#x = solve(a,b)


def f1(x):
    return 1


def f2(x):
    return x**2


def solveBVP1(f,N):
    h = 1.0/N
    x = zeros(N-1)
    for i in range(N-1):
        x[i] = (i+1)*h
    f_h = [f(num) for num in x]
    
    diag_h = [2]*(N-1)
    sub_diag_h = [-1]*(N-2)
    a_h = array([[0]+sub_diag_h, diag_h])
    
    u_h = solveh_banded(a_h,f_h)
    figure()
    plot(x,u_h)
    xlabel('x')
    ylabel('u')
    title('n='+ str(N))
    show()
    
    

def solveBVP2(f,N):
    h = 1.0/N
    x = zeros(N-1)
    for i in range(N-1):
        x[i] = (i+1)*h
    f_h = [f(num) for num in x]
    
    diag_h = [2+h**2]*(N-1)
    sub_diag_h = [-1]*(N-2)
    a_h = array([[0]+sub_diag_h, diag_h])
    
    u_h = solveh_banded(a_h,f_h)
    figure()
    plot(x,u_h)
    xlabel('x')
    ylabel('u')
    title('n='+ str(N))
    show()








