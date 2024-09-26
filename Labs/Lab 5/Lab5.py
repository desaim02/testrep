import numpy as np
import matplotlib.pyplot as plt


# 1

#All x_0 such that the absolute value of g'(x_0) where g' is defined by Newton's Method is said to be in the basin of convergence

# 2 and 4 (made changes to same function defintion)
def bisection(f,a,b,tol):

#     first verify there is a root we can find in the interval 
    count = 0
    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,count]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,count]
    
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier, count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier, count]

def newton(f,fp,p0,tol,Nmax):
  
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1)
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def newt_bisection(f,fp,fpp,a,b,tol):

#     first verify there is a root we can find in the interval 
    count = 0
    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,count]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,count]
    

    d = 0.5*(a+b) #d is midpoint 
    gp = lambda x: (f(x)*fpp(x)) / (fp(x))**2
    if abs(gp(d)) <1:
        print ('midpoint is in the basin of convergence')
        Nmax = 100
        [p,pstar,info,it] = newton(f,fp,d,tol,Nmax)
        print('the array of the iterates is', p)
        print('the last iterate is', pstar)
        print('the sucess message (0 good) is', info)
        print('the number of iterations is', it)
        return [pstar]
    
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier, count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier, count]

f = lambda x: np.exp(x**2+7*x-30) - 1
fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)
fpp = lambda x: (2*x+7)*(2*x+7)*np.exp(x**2+7*x-30) + 2*np.exp(x**2+7*x-30)

def driver6a():
   f = lambda x: np.exp(x**2+7*x-30) - 1
   a = 2
   b = 4.5
   tol = 10**-8
   [astar, ier,count] = bisection(f,a,b,tol)
   print('the astar is', astar)
   print('the ier is', ier)
   print('the count is', count)

driver6a()

def driver6b():
    f = lambda x: np.exp(x**2+7*x-30) - 1
    fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)
    p0 = 4.5
    tol = 10**-8
    Nmax = 100
    [p,pstar,info,it] = newton(f,fp,p0,tol,Nmax)
    print('the array of the iterates is', p)
    print('the last iterate is', pstar)
    print('the sucess message (0 good) is', info)
    print('the number of iterations is', it)

driver6b()

def driver6c():
    a = 2
    b = 4.5
    tol = 10**-8
    f = lambda x: np.exp(x**2+7*x-30) - 1
    fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)
    fpp = lambda x: (2*x+7)*(2*x+7)*np.exp(x**2+7*x-30) + 2*np.exp(x**2+7*x-30)
    newt_bisection(f,fp,fpp,a,b,tol)

driver6c()

#The bisection method took 27 iterations. Newton's method took 26 iterations, and the hybrid method took 7 iterations. Our hybrid method converged signicantly faster than the rest. 

   
