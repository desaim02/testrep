#old fixed point code

# import libraries
import numpy as np

def modfixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    x = np.array(x0)
    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       x = np.append(x,x1)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          print(x)
          return [xstar,x, ier]
       x0 = x1

    xstar = x1
    ier = 1
    print(x)
    return [xstar, x, ier]
    
# 2.2 Exercises 

# 1

def compute_order(x,xstar):

    diff1 = np.abs(x[1::] -xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)

    _lambda = np.exp(fit[1])
    alpha = fit[0]
    print ('lambda is', _lambda)
    print('alphda is', alpha)
    return fit

#2
#a

def driver2a():
     f1 = lambda x: (10/(x+4))**.5
     Nmax = 100
     tol = 1e-10
     x0 = 1.5
     [xstar,x, ier] = modfixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('x is', x)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
     print('the amount of iterations it takes for the FPI to converge with the abosulte tol is ', len(x) )

driver2a()

#b
def driver2b():
     f1 = lambda x: (10/(x+4))**.5
     Nmax = 100
     tol = 1e-10
     x0 = 1.5
     [xstar,x,ier] = modfixedpt(f1,x0,tol,Nmax)
     xstar = x[-1]
     x = x-x[-1]
     compute_order(x,xstar)

driver2b()

#3
# def Aitken(x):
#     pn = np.abs(x[0::])
#     pn1 = np.abs(x[1::])
#     pn2 = np.abs(x[2::])
#     p = pn- ((pn1 - pn)**2 / pn2-2*pn1+pn)
#     print(p)

def Aitken(x):
    xn = x[:-2]
    xn1 = x[1:-1]
    xn2 = x[2:]
    return xn- ((xn1 - xn)**2 / xn2-2*xn1-xn)

def driver3a():
     f1 = lambda x: (10/(x+4))**.5
     Nmax = 100
     tol = 1e-10
     x0 = 1.5
     [xstar,x,ier] = modfixedpt(f1,x0,tol,Nmax)
     xstar = x[-1]
     x = x-x[-1]
     aitkensx = Aitken(x)
     compute_order(aitkensx, xstar)
    
driver3a()
