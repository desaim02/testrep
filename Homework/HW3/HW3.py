# import libraries
import numpy as np
import matplotlib.pyplot as plt

#functions
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

def fixedpt(f,x0,tol,Nmax):
    # x0 = initial guess 
    # Nmax = max number of iterations'''
    # tol = stopping tolerance'''
    count = 0
    while (count <Nmax):
        count = count +1
        x1 = f(x0)
        if (abs(x1-x0) <tol):
            xstar = x1
            ier = 0
            return [xstar,ier]
        x0 = x1
    xstar = x1
    ier = 1
    return [xstar, ier]

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
    


#Exercise 1

#show plot so we can intuitively see where root is / find a range for it
f = lambda x: 2*x-1-np.sin(x)
xvals = np.linspace(-1,2,100)
yvals = f(xvals)

plt.plot(xvals,yvals)
plt.title('2x-1-sin(x)')
plt.axhline(y=0)
plt.axvline(x=0, linestyle='--')
plt.axvline(x=np.pi/2,linestyle='--')
plt.show()

#1c
def driver1c():
# use routines    
    f = lambda x:  2*x-1-np.sin(x)
    tol = 10e-8
    a = 0
    b = 1

    [astar,ier,count] = bisection(f,a,b,tol)

    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('the number of iterations was' , count)

driver1c()

#Exercise 2
#2

def driver2a():
    f = lambda x: (x-5)**9
    tol = 1e-4
    a = 4.82
    b = 5.2

    [astar,ier,count] = bisection(f,a,b,tol)

    print('for part a:')
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('the number of iterations was' , count)

driver2a()

def driver2b():
    f = lambda x: x**9-45*x**8+900*x**7-10500*x**6+78750*x**5-393750*x**4+1312500*x**3-2812500*x**2+3515625*x-1953125
    tol = 1e-4
    a = 4.82
    b = 5.2

    [astar,ier,count] = bisection(f,a,b,tol)

    print('for part b:')
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('the number of iterations was' , count)

driver2b()

f1 = lambda x: (x-5)**9
f2 = lambda x: x**9-45*x**8+900*x**7-10500*x**6+78750*x**5-393750*x**4+1312500*x**3-2812500*x**2+3515625*x-1953125
xvals = np.linspace(4.8,5.22,100)
yvals1 = f1(xvals)
yvals2 = f2(xvals)

plt.plot(xvals,yvals1, color='orange')
plt.plot(xvals,yvals2, color='green')
plt.axhline(y=0, linestyle='--')
plt.show()

#Exercise 3
def driver3b():
   
  f = lambda x: x**3+x-4
  tol = 10e-3
  a = 1
  b = 4

  [astar,ier,count] = bisection(f,a,b,tol)

  print('the approximate root is',astar)
  print('the error message reads:',ier)
  print('f(astar) =', f(astar))
  print('the number of iterations was' , count)

  driver3b()

#Exercise 4

def compute_order(x,xstar):

    diff1 = np.abs(x[1::] -xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)

    _lambda = np.exp(fit[1])
    alpha = fit[0]
    print ('lambda is', _lambda)
    print('alphda is', alpha)
    return fit

def driver4c():
  f = lambda x: 12/(1+x)
  x0 = 2.9
  tol =10e-10
  nmax=100
  [xstar,x,ier] = modfixedpt(f,x0,tol,nmax)
  xstar = x[-1]
  x = x-x[-1]

  compute_order(x,xstar)

driver4c()




#Exercise 5

def driver5a():
  f = lambda x: x-4*np.sin(2*x)-3
  xvals = np.linspace(-1,2*np.pi,100)
  yvals= f(xvals)

  plt.plot(xvals,yvals)
  plt.title('x-4sin(2x)-3=0')
  plt.axhline(y=0)
  plt.show()

driver5a()

def driver5b():
  f = lambda x: -np.sin(2*x)+5*(x/4)-3/4
  g = lambda x: -2*np.cos(2*x) + (5/4)
  tol = .5*10**-10 #because we want 10 correct digits 
  nmax = 100
  x01 = -.9
  x02 = -.5
  x03 = 1.4
  x04 = 3
  x05 = 4.5


  [xstar1,ier1] = fixedpt(f,x01,tol,nmax)
  print ('xstar is ', xstar1)
  print ('g is', g(x01))
  print('ier is', ier1)

  [xstar2,ier2] = fixedpt(f,x02,tol,nmax)
  print ('xstar is ', xstar2)
  print ('g is', g(x02))
  print('ier is', ier2)

  [xstar3,ier3] = fixedpt(f,x03,tol,nmax)
  print ('xstar is ', xstar3)
  print ('g is', g(x03))
  print('ier is', ier3)


  [xstar4,ier4] = fixedpt(f,x04,tol,nmax)
  print ('xstar is ', xstar4)
  print ('g is', g(x04))
  print('ier is', ier4)


  [xstar5,ier5] = fixedpt(f,x05,tol,nmax)
  print ('xstar is ', xstar5)
  print ('g is', g(x05))
  print('ier is', ier5)

driver5b()









