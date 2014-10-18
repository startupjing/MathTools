from __future__ import division
from pylab import *

def newtonFractal(f, Jf, n=32, p=512):
    """Display a Newton fractal for the function f with Jacobian Jf.
    n is the number of Newton steps (default: 32), and
    p is the width and height of the image in pixels (default: 512)."""

    # Create a p-by-p grid of points in the square [-1,1] x [-1,1].
    x = array(meshgrid(linspace(-1,1,p), linspace(-1,1,p)))

    # Do n Newton steps starting from each point in the grid.
    for k in range(n):
        # The following line is just a clever way of doing a Newton
        # step on the entire array of grid points simultaneously.
        # Since NumPy array operations are highly optimized, this is
        # much faster than doing a "for loop" over the grid points.
        x = x - solve(transpose(Jf(x),(3,2,0,1)), f(x).T).T

    # Show a p-by-p pixel image, where each starting point is colored
    # according to the angle of the point in the plane where it ended
    # up. This angle is given by arctan2(*x) = arctan(x1/x0).
    imshow(arctan2(*x), extent=[-1,1,-1,1])
    
def f(x):
    return array([x[0]**3 - 3*x[0]*x[1]**2 - 1, 3*x[0]**2*x[1] - x[1]**3])

def Jf(x):
    return array([[3*x[0]**2-3*x[1]**2, -6*x[0]*x[1]], [6*x[0]*x[1], 3*x[0]**2-3*x[1]**2]])    

def g(x):
    return array([x[0]**4+x[1]**4-6*(x[0]**2)*(x[1]**2)-1, 4*(x[0]**3)*x[1]-4*x[0]*(x[1]**3)])

def Jg(x):
    return array([[4*(x[0]**3)-12*x[0]*(x[1]**2), 4*(x[1]**3)-12*(x[0]**2)*x[1]],[-4*(x[1]**3)+12*(x[0]**2)*x[1], 4*(x[0]**3)-12*x[0]*(x[1]**2)]])

def hilbert(n):
    """create a n by n Hilbert matrix"""
    hilbert = numpy.zeros(shape=(n,n))
    for row in range(n):
        for col in range(n):
            hilbert[row, col] = 1/(row+col+1);
    return hilbert
    
def compare():
    """solve problem 3(b)"""
    a = hilbert(20);
    """x = numpy.zeros(shape=(20,1));"""
    x = zeros(20);
    for i in range(20):
        x[i] = i+1
    b = dot(a,x);
    x_soln = solve(a,b);
    delta_x = x - x_soln;
    delta_b = dot(a,x+delta_x) - b;
    ratio1 = norm(delta_b)/norm(b);
    ratio2 = norm(delta_x)/norm(x);
    print 'difference_b: '
    print ratio1
    print 'difference_x: '
    print ratio2
    
    

    
        
            
