# import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.special as sp
import pandas as pd



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

def compute_order(x,xstar):

    diff1 = np.abs(x[1::] -xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)

    _lambda = np.exp(fit[1])
    alpha = fit[0]
    print ('lambda is', _lambda)
    print('alphda is', alpha)
    return fit

def secant(f,x0,x1,tol,Nmax):
    """
  Secant iteration.
  
  Inputs:
    f - function 
    x0   - initial guess 1
    x1 - inital guess 2
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
    p[0] = x0
    p[1] = x1
    for it in range(Nmax):
        x2 =  x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        p[it+2] = x2
        if (abs(x2-x1) < tol):
            pstar = x2
            info = 0
            return [p,pstar,info,it]
        x0 = x1
        x1=x2
    pstar = x2
    info = 1
    return [p,pstar,info,it]

#1

#a
T_i = 20
T_s = -15
alpha = .138*10**-6

def driver1a():
    x_bar = 2
    t = 5184000
    f = lambda x: sp.erf(x / (2 * np.sqrt(alpha * t))) * (T_i - T_s) + T_s


    print('f(2) = ', f(x_bar) )
    xvals = np.linspace(0,x_bar,100)
    yvals = f(xvals)
    plt.plot(xvals,yvals)
    plt.axhline(y=0, linestyle='--')
    plt.show()

driver1a()

def driver1b():
   t = 5184000
   f = lambda x: sp.erf(x / (2 * np.sqrt(alpha * t))) * (T_i - T_s) + T_s
   a_0 = 0
   b_0 = 2
   tol = 10**-13
   [astar,ier,count] = bisection(f,a_0,b_0,tol)
   print('the astar is', astar)
   print('the ier is', ier)
   print('the count is', count)

driver1b()

def driver1c():
   t = 5184000
   f = lambda x: sp.erf(x / (2 * np.sqrt(alpha * t))) * (T_i - T_s) + T_s
   fp = lambda x: ((T_i - T_s) / np.sqrt(np.pi * alpha * t)) * np.exp(- (x / (2 * np.sqrt(alpha * t)))**2)

   Nmax = 10 #can adjust if neccessary
   tol = 10**-13
   p0 = .01
   p01 = 2
   
   [p,pstar,info,it] = newton(f,fp,p0,tol,Nmax)
   print ("For x0 =.01:")
   print('the array of the iterates is', p)
   print('the last iterate is', pstar)
   print('the sucess message (0 good) is', info)
   print('the number of iterations is', it)

   [p,pstar,info,it] = newton(f,fp,p01,tol,Nmax)
   print ("For x0 =xbar:")
   print('the array of the iterates is', p)
   print('the last iterate is', pstar)
   print('the sucess message (0 good) is', info)
   print('the number of iterations is', it)

driver1c()

#4

def driver4i():
    print ('the Newtons method is')
    f = lambda x: np.exp(3*x)-27*x**6+27*x**4*np.exp(x)-9*x**2*np.exp(2*x)
    fp  = lambda x:3*np.exp(3*x)-162*x**5+27*x**4*np.exp(x)+108*x**3*np.exp(x)-18*x**2*np.exp(2*x)-18*x*np.exp(2*x)
    tol = 10**-8
    Nmax=35
    p0=4
    [p,pstar,info,it] = newton(f,fp,p0,tol,Nmax)
    print('the array of the iterates is', p)
    print('the last iterate is', pstar)
    print('the sucess message (0 good) is', info)
    print('the number of iterations is', it)

    p = p[:-5] #editing iteration array so that it doesn't contain 0s or identified root 
    compute_order(p,pstar)

driver4i()


def driver4ii():
    print ('the modified Newtons method from class is')
    f = lambda x: np.exp(3*x)-27*x**6+27*x**4*np.exp(x)-9*x**2*np.exp(2*x)
    fp  = lambda x:3*np.exp(3*x)-162*x**5+27*x**4*np.exp(x)+108*x**3*np.exp(x)-18*x**2*np.exp(2*x)-18*x*np.exp(2*x)
    fpp = lambda x: 9*np.exp(3*x)-810*x**4+27*x**4*np.exp(x)+216*x**3*np.exp(x)+324*x**2*np.exp(x)-36*x**2*np.exp(2*x)-72*x*np.exp(2*x)-18*np.exp(2*x)
    g = lambda x: x - (f(x)/fp(x))
    gp = lambda x: (f(x)*fpp(x)) / (fp(x))**2
    tol = 10**-8
    Nmax=100
    p0=4
    [p,pstar,info,it] = newton(g,gp,p0,tol,Nmax)
    print('the array of the iterates is', p)
    print('the last iterate is', pstar)
    print('the sucess message (0 good) is', info)
    print('the number of iterations is', it)

    p = p[:-1]
    compute_order(p,pstar)

driver4ii()

def driver4iii():
  print ('the modified Newtons method from part 2c is')
  f = lambda x: np.exp(3*x)-27*x**6+27*x**4*np.exp(x)-9*x**2*np.exp(2*x)
  fp  = lambda x:3*np.exp(3*x)-162*x**5+27*x**4*np.exp(x)+108*x**3*np.exp(x)-18*x**2*np.exp(2*x)-18*x*np.exp(2*x)
  fpp = lambda x: 9*np.exp(3*x)-810*x**4+27*x**4*np.exp(x)+216*x**3*np.exp(x)+324*x**2*np.exp(x)-36*x**2*np.exp(2*x)-72*x*np.exp(2*x)-18*np.exp(2*x)
  tol = 10**-8
  Nmax=100
  p0=4
   #code below was to decide which value to use for m 
  # for i in range(1,5): #trying out this Newton's method with varing m's to see which is the most successful 
  #   m = i
  #   g = lambda x: x-m*(f(x)/fp(x))
  #   [xstar,x,ier] = modfixedpt(g,p0,tol,Nmax)
  #   print (" m = ", i)
  #   print('the array of the iterates is', x)
  #   print('the last iterate is', xstar)
  #   print ('the number of iteration is', len(x))
  #   print('the success message (0 good) is', ier)

    #Note: m=3 was the first one that converged and also did so with less iterations than m=4. 

  m=3
  g = lambda x: x-m*(f(x)/fp(x))
  [xstar,x,ier] = modfixedpt(g,p0,tol,Nmax)
  print('the array of the iterates is', x)
  print('the last iterate is', xstar)
  print ('the number of iteration is', len(x))
  print('the success message (0 good) is', ier)
  x = x[:-2]
  print(x)
  compute_order(x,xstar)

driver4iii()


#5


def drivernewt5():
  f=lambda x: x**6-x-1
  fp = lambda x: 6*x**5-1
  p0 = 2
  tol = 10**-8
  Nmax=10
  [p,pstar,info,it] = newton(f,fp,p0,tol,Nmax)
  print('the array of the iterates is', p)
  print('the last iterate is', pstar)
  print('the sucess message (0 good) is', info)
  print('the number of iterations is', it)

  #Part A: Create Error Table 
  it_num= range(len(p))
  iterate = [p[i] for i in it_num]
  it_err = [p[i]-pstar for i in it_num]

  df = pd.DataFrame({
    'Iteration': list(it_num),
    'Iteration Guess': iterate,
    'Error': it_err
    })
  print("Newton's Method Error Table:")
  print(df)
  

  # Part B

  err = abs(p-pstar)
  print(err)

  plt.plot(range(len(err)), np.log(err))
  plt.title('Newton Method')
  plt.xlabel("Iteration")
  plt.ylabel("Log of Error")
  plt.show()  

  newt_x = range(len(err))
  newt_y = np.log(err)
  return [newt_x,newt_y]
  # it_err_array = np.array(it_err)
  # it_err_array =  it_err_array[:-3]
  
  # xvals1 = np.abs(it_err_array[1:])
  # xvals2 = np.abs(it_err_array[:-1])
  # print(xvals2)

  # plt.loglog(xvals1,xvals2, color='green')
  # # plt.xlabel("|x_{k+1} - alpha|")
  # # plt.ylabel("|x_k - alpha|") 
  # # plt.legend()
  # plt.show()

#drivernewt5()

def driversec5():
    f=lambda x: x**6-x-1
    p0 = 2
    p1 = 1
    tol = 10**-8
    Nmax=10
    [p,pstar,info,it] = secant(f,p0,p1,tol,Nmax)
    print('the array of the iterates is', p)
    print('the last iterate is', pstar)
    print('the sucess message (0 good) is', info)
    print('the number of iterations is', it)

    #Part A: Create Error Table 
    it_num= range(len(p))
    iterate = [p[i] for i in it_num]
    it_err = [p[i]-pstar for i in it_num]

    df = pd.DataFrame({
      'Iteration': list(it_num),
      'Iteration Guess': iterate,
      'Error': it_err
      })
    print("Secant Method Error Table:")
    print(df)


    #Part B: Plot Error
    err = abs(p-pstar)
    print(err)

    plt.plot(range(len(err)), np.log(err))
    plt.title('Secant Method')
    plt.xlabel("Iteration")
    plt.ylabel("Log of Error")
    plt.show()

    sec_x = range(len(err))
    sec_y = np.log(err)
    return [sec_x,sec_y]

#driversec5()


def driver5():
    [sec_x,sec_y] = driversec5()
    [newt_x,newt_y] = drivernewt5()


    plt.plot(newt_x, newt_y, color='green')
    plt.plot(sec_x,sec_y,color='orange')
    plt.title('Newton vs. Secant Method')
    plt.xlabel("Iteration")
    plt.ylabel("Log of Error")
    plt.legend()
    plt.show()

driver5()

   