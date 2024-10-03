import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def compute_order(x,xstar):

    diff1 = np.abs(x[1::] -xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)

    _lambda = np.exp(fit[1])
    alpha = fit[0]
    print ('lambda is', _lambda)
    print('alphda is', alpha)
    return fit


#Pre-Lab

s = np.pi/2
f = lambda x: np.cos(x)
h = .01 * 2.0 ** (-np.arange(0,10))
tol = 1**-8

def forward(s,f,h):
    return (f(s+h)-f(s))/h

def centered(s,f,h):
    return (f(s+h)-f(s-h))/(2*h)

print("The result for the forward difference is: ")
print(forward(s,f,h))
print("The order of convergence is: ")
compute_order(forward(s,f,h)[:-1],forward(s,f,h)[-1])

print("The result for the centered difference is: ")
print(centered(s,f,h))
print("The order of convergence is: ")
compute_order(centered(s,f,h)[:-1],centered(s,f,h)[-1])


def LazyNewton(evalF, evalJ, x0,tol,Nmax=500):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
        if (its == 0 or its%5 == 0):
            J = evalJ(x0)
            Jinv = inv(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
        
        x0 = x1
    
    xstar = x1
    ier = 1 
    print( "ran out of iterations")
    return[xstar,ier,its]

# Lab 

# 3.2: Slacker Newton 

# 1: I'm going to recompute the Jacobian every 5 iterates. 

# 2: Where is lazy newton code

def drivera():

    def evalF(x): #function 
        F = np.zeros(2)
        F[0] = 4*(x[0]**2)+(x[1]**2)-4
        F[1] = x[0] + x[1] - np.sin(x[0]-x[1])
        return F

    def evalJ(x): #hve to change for jacobian
        J = np.array([[8 *x[0], 2*x[1]],
            [1-np.cos(x[0]-x[1]), 1+np.cos(x[0]-x[1])]])
        return J
    
    x0 = np.array([1,0])
    tol = 10**-10
    
    [xstar, ier,its] = LazyNewton(evalF, evalJ, x0,tol)

    print('Lazy Newton: the error message reads:',ier)
    print('xstar is:',xstar)
    print('Lazy Newton: number of iterations is:',its)


drivera()


# My iteration, which computes the Jacobian every 5 iterations, converges in 6 total iterations. 

