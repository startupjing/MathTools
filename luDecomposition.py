from __future__ import division
from pylab import *


def lu(A):
    """LU decomposition of A without pivoting (Doolittle's method)."""
    n = shape(A)[0]
    L = eye(n)
    U = zeros([n,n])
    for k in range(n):
        U[k,k:] = A[k,k:] - dot(L[k,:k], U[:k,k:])        
        L[k+1:,k] = (A[k+1:,k] - dot(L[k+1:,:k], U[:k,k]))/U[k,k]
    return L, U

def uSolve(U,b):
    """Solve an upper triangular system Ux = b by back substitution."""
    n = size(b)
    x = zeros(n)
    # the reversed() function lets you iterate backwards over a list
    for i in reversed(range(n)):
        x[i] = (b[i] - dot(U[i,i+1:], x[i+1:]))/U[i,i]
    return x
    
def lSolve(L,b):
    """Solve a lower triangular system Lx = b by forward substitution."""
    n = size(b)
    x = zeros(n)
    for i in range(n):
        x[i] = (b[i] - dot(L[i,0:i], x[0:i]))/L[i,i]
    return x

def luSolve(L,U,b):
    """Solve system LUb by forward and backward substituion.
       Solve L(Ux) = b first to get Ux = temp using lsolve(),
       then solve Ux = temp using usolve() to get x"""
    temp = lSolve(L,b)
    return uSolve(U,temp)
    
def explore(n):
    """Find smallest k such that error occurs
       using the range(n), print message if a larger
       range is needed"""
    for k in range(n):
      test = array([[math.pow(10,-k),1],[1,1]])
      L,U = lu(test)
      result = dot(L,U)
      if(result[1][1] == 0):
          return k
    return "need larger range"

def verify(k):
    """verify no error occurs with pivot matrix using k"""
    pivotA = array([[1,1],[math.pow(10,-k),1]])
    L,U = lu(pivotA)
    return pivotA, dot(L,U)
    
    
