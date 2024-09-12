# import libraries
import numpy as np
import matplotlib.pyplot as plt

def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
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
    return [astar, ier]

#Exercise 1


def driver1():
# use routines    
    f = lambda x: x**2*(x-1)
    tol = 10e-5

    [astara,iera] = bisection(f,.5,2,tol)

    print('for part a:')
    print('the approximate root is',astara)
    print('the error message reads:',iera)
    print('f(astar) =', f(astara))

    [astarb,ierb] = bisection(f,-1,.5,tol)

    print('for part b:')
    print('the approximate root is',astarb)
    print('the error message reads:',ierb)
    print('f(astar) =', f(astarb))

    [astarc,ierc] = bisection(f,-1,2,tol)

    print('for part c:')
    print('the approximate root is',astarc)
    print('the error message reads:',ierc)
    print('f(astar) =', f(astarc))

driver1()

#show plot so we can intuitively see where roots are
f = lambda x: x**2*(x-1)
xvals = np.linspace(-1,2,100)
yvals = f(xvals)

plt.plot(xvals,yvals)
plt.axhline(y=0)
plt.show()

# Question 1 Analysis/Response
# The bisection method is successful for parts a and c, but not b. This is because parts a and c successfully fufull the bisection method requirements of having a postive and negative point. Whereas for part b, even though there is a visual roots within (-1,.5), the points left and right of the root are both negative. Thus, the bisection method didn't work to identify the root. 


def driver2a():
    f = lambda x: (x-1)*(x-3)*(x-5)
    tol = 10e-5
    a = 0 
    b = 2.4

    [astar,ier] = bisection(f,a,b,tol)

    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))


    xvals =  np.linspace(a,b,100)
    yvals = f(xvals)

    plt.plot(xvals,yvals)
    plt.axhline(y=0)
    plt.show()  

driver2a()

# Question 2a
#The code was succesful for this question as points a and b and the function fulfilled the requirements (f is continous over the interval (a,b) and a*b = negative (i.e. they are opposite signs))

def driver2b():
    f = lambda x: (x-1)**2 * (x-3)
    tol = 10e-5
    a = 0 
    b = 2

    [astar,ier] = bisection(f,a,b,tol)

    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))


    xvals =  np.linspace(a,b,100)
    yvals = f(xvals)

    plt.plot(xvals,yvals)
    plt.axhline(y=0)
    plt.show()  

driver2b()

#Question 2b 
#This code was unsuccessful because our points (0,2) were both negative, and hence did not fulfull our bisection method requirements. 

def driver2c1():
    f = lambda x: np.sin(x)
    tol = 10e-5
    a = 0 
    b = .1

    [astar,ier] = bisection(f,a,b,tol)

    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))


    xvals =  np.linspace(a,b,100)
    yvals = f(xvals)

    plt.plot(xvals,yvals)
    plt.axhline(y=0)
    plt.show()  

driver2c1()

# Question 2c
#For the first prt of 2c with the inteval (0,.1), it just so happened that our code worked.. However, our interval did not fulfull the bisection method requirements, so we should not have used it in this scenario. 

def driver2c2():
    f = lambda x: np.sin(x)
    tol = 10e-5
    a = .5
    b = (3*np.pi)/4

    [astar,ier] = bisection(f,a,b,tol)

    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))


    xvals =  np.linspace(a,b,100)
    yvals = f(xvals)

    plt.plot(xvals,yvals)
    plt.axhline(y=0)
    plt.show()  

driver2c2()



# For the second part of 2c, our code did not work. This is because we, once again, did not have opposite signs. In this case, there is no point in between our interval (.5, 3pi.4) that is a root. 

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

# Question 3 
# Recall: A pt,c, is a fixed point of f(x) if f(c) = c. 

f3a = lambda x: x * (1+(7-(x**5))/(x**2))**3
f3b = lambda x: x - (x**5-7)/x**2
f3c = lambda x: x - (x**5-7)/(5*x**4)
f3d = lambda x: x - (((x**5)-7)/12)

x = 7**(1/5)

print(f3a(x))
print(f3b(x))
print(f3c(x))
print(f3d(x))
print("x= ", 7**(1/5))

Nmax = 100
tol = 1e-10
x0 = 1

[xstar,ier] = fixedpt(f3a,x0,tol,Nmax)
print('part a :')
print('the approximate fixed point is:',xstar)
print('f3a(xstar):',f3a(xstar))
print('Error message reads:',ier)

print('part b :')
[xstar,ier] = fixedpt(f3b,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f3b(xstar):',f3b(xstar))
print('Error message reads:',ier)

print('part c :')
[xstar,ier] = fixedpt(f3c,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f3c(xstar):',f3c(xstar))
print('Error message reads:',ier)

print('part d :')
[xstar,ier] = fixedpt(f3d,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f3d(xstar):',f3d(xstar))
print('Error message reads:',ier)

#It converged for part c and d. It doesn't converge for parts a and b because the slope is greater than 1. 