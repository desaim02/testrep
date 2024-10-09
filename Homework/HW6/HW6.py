import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]
           
def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]   

def Broyden(x0,tol,Nmax):
    '''tol = desired accuracy
    Nmax = max number of iterations'''

    '''Sherman-Morrison 
   (A+xy^T)^{-1} = A^{-1}-1/p*(A^{-1}xy^TA^{-1})
    where p = 1+y^TA^{-1}Ax'''

    '''In Newton
    x_k+1 = xk -(G(x_k))^{-1}*F(x_k)'''

    '''In Broyden 
    x = [F(xk)-F(xk-1)-\hat{G}_k-1(xk-xk-1)
    y = x_k-x_k-1/||x_k-x_k-1||^2'''

    ''' implemented as in equation (10.16) on page 650 of text'''
    
    '''initialize with 1 newton step'''
    
    A0 = evalJ(x0)

    v = evalF(x0)
    A = np.linalg.inv(A0)

    s = -A.dot(v)
    xk = x0+s
    for  its in range(Nmax):
       '''(save v from previous step)'''
       w = v
       ''' create new v'''
       v = evalF(xk)
       '''y_k = F(xk)-F(xk-1)'''
       y = v-w;                   
       '''-A_{k-1}^{-1}y_k'''
       z = -A.dot(y)
       ''' p = s_k^tA_{k-1}^{-1}y_k'''
       p = -np.dot(s,z)                 
       u = np.dot(s,A) 
       ''' A = A_k^{-1} via Morrison formula'''
       tmp = s+z
       tmp2 = np.outer(tmp,u)
       A = A+1./p*tmp2
       ''' -A_k^{-1}F(x_k)'''
       s = -A.dot(v)
       xk = xk+s
       if (norm(s)<tol):
          alpha = xk
          ier = 0
          return[alpha,ier,its]
    alpha = xk
    ier = 1
    return[alpha,ier,its]


#note: x[0] = x and x[1] = y
def evalF(x): #function 
        F = np.zeros(2)
        F[0] = x[0]**2 + x[1]**2 -4
        F[1] = np.exp(x[0]) +x[1] -1
        return F

#  jacobian is 2x, 2y
# e^x, 1
def evalJ(x): #hve to change for jacobian
    J = np.array([[2 *x[0], 2*x[1]],
        [np.exp(x[0]), 1]])
    return J


def driverBroyden():
    print("The Broyden Results are")
    x1 = np.array([1,1])
    x2 = np.array([1,-1])
    x3 = np.array([0,0])
    tol = 10**-10
    Nmax =100

    t = time.time()
    for j in range(20):
      [xstar1,ier1,its1] = Broyden(x1, tol,Nmax)     
    elapsed1 = time.time()-t
    print(xstar1)
    print('Broyden: the error message reads:',ier1)
    print('Broyden: took this many seconds:',elapsed1/20)
    print('Broyden: number of iterations is:',its1)

    t = time.time()
    for j in range(20):
      [xstar2,ier2,its2] = Broyden(x2, tol,Nmax)     
    elapsed2 = time.time()-t
    print(xstar2)
    print('Broyden: the error message reads:',ier2)
    print('Broyden: took this many seconds:',elapsed2/20)
    print('Broyden: number of iterations is:',its2)


    # t = time.time()
    # for j in range(20):
    #   [xstar3,ier3,its3] = Broyden(x3, tol,Nmax)     
    # elapsed3 = time.time()-t
    # print(xstar3)
    # print('Broyden: the error message reads:',ier3)
    # print('Broyden: took this many seconds:',elapsed3/20)
    # print('Broyden: number of iterations is:',its3)

    # print ("The results are")
    # return (xstar1,xstar2)
    #return (xstar1,xstar2,xstar3)

    #print("The iterations are")
    #return (its1,its2,its3)

    #print("The timing is")
    #return (elapsed1,elapsed2,elpased3)

def driverLazy():
    print("Lazy Newton Results:")
    x1 = np.array([1,1])
    x2 = np.array([1,-1])
    x3 = np.array([0,0])
    tol = 10**-10
    Nmax =100

    t = time.time()
    for j in range(20):
      [xstar1,ier1,its1] =  LazyNewton(x1,tol,Nmax)
    elapsed1 = time.time()-t
    print(xstar1)
    print('Lazy Newton: the error message reads:',ier1)
    print('Lazy Newton: took this many seconds:',elapsed1/20)
    print('Lazy Newton: number of iterations is:',its1)

    t = time.time()
    for j in range(20):
      [xstar2,ier2,its2] =  LazyNewton(x2,tol,Nmax)
    elapsed2 = time.time()-t
    print(xstar2)
    print('Lazy Newton: the error message reads:',ier2)
    print('Lazy Newton: took this many seconds:',elapsed2/20)
    print('Lazy Newton: number of iterations is:',its2)
    
    # t = time.time()
    # for j in range(20):
    #   [xstar3,ier3,its3] =  LazyNewton(x3,tol,Nmax)
    # elapsed3 = time.time()-t
    # print(xstar3)
    # print('Lazy Newton: the error message reads:',ier3)
    # print('Lazy Newton: took this many seconds:',elapsed3/20)
    # print('Lazy Newton: number of iterations is:',its3)

    # print ("The results are")
    #return (xstar1,xstar2)
    #return (xstar1,xstar2,xstar3)

    # print("The iterations are")
    # return (its1,its2)

    # print("The timing is")
    # return (elapsed1,elapsed2)


def driverNewton():
    print("The Newton Results are")
    x1 = np.array([1,1])
    x2 = np.array([1,-1])
    x3 = np.array([0,0])
    tol = 10**-10
    Nmax =100

    t = time.time()
    for j in range(50):
      [xstar1,ier1,its1] =  Newton(x1,tol,Nmax)
    elapsed1 = time.time()-t
    print(xstar1)
    print('Newton: the error message reads:',ier1) 
    print('Newton: took this many seconds:',elapsed1/50)
    print('Netwon: number of iterations is:',its1)

    t = time.time()
    for j in range(50):
      [xstar2,ier2,its2] =  Newton(x2,tol,Nmax)
    elapsed2 = time.time()-t
    print(xstar2)
    print('Newton: the error message reads:',ier2) 
    print('Newton: took this many seconds:',elapsed2/50)
    print('Netwon: number of iterations is:',its2)

    # t = time.time()
    # for j in range(50):
    #   [xstar3,ier3,its3] =  Newton(x3,tol,Nmax)
    # elapsed3 = time.time()-t
    # print(xstar3)
    # print('Newton: the error message reads:',ier3) 
    # print('Newton: took this many seconds:',elapsed3/50)
    # print('Netwon: number of iterations is:',its3)

    # print ("The results are")
    # return (xstar1,xstar2)
    #return (xstar1,xstar2,xstar3)

    #print("The iterations are")
    #return (its1,its2,its3)

    #print("The timing is")
    #return (elapsed1,elapsed2,elpased3)


def driver():

    driverNewton()
    driverLazy()
    driverBroyden()

driver()

# Question 2 
tol =  10 ** -6
Nmax = 100
x = [.5,.5,.5]


#functions:
def evalF(x):

    F = np.zeros(3)
    F[0] = x[0] +np.cos(x[0]*x[1]*x[2])-1
    F[1] = (1-x[0])**(0.25) + x[1] +0.05*x[2]**2 -0.15*x[2]-1
    F[2] = -x[0]**2-0.1*x[1]**2 +0.01*x[1]+x[2] -1
    return F

def evalJ(x): 

    J =np.array([[1+x[1]*x[2]*np.sin(x[0]*x[1]*x[2]),x[0]*x[2]*np.sin(x[0]*x[1]*x[2]),x[1]*x[0]*np.sin(x[0]*x[1]*x[2])],
          [-0.25*(1-x[0])**(-0.75),1,0.1*x[2]-0.15],
          [-2*x[0],-0.2*x[1]+0.01,1]])
    return J

def evalg(x):

    F = evalF(x)
    g = F[0]**2 + F[1]**2 + F[2]**2
    return g

def eval_gradg(x):
    F = evalF(x)
    J = evalJ(x)
    
    gradg = np.transpose(J).dot(F)
    return gradg


###############################
### steepest descent code

def SteepestDescent(x,tol,Nmax):
    
    for its in range(Nmax):
        g1 = evalg(x)
        z = eval_gradg(x)
        z0 = norm(z)

        if z0 == 0:
            print("zero gradient")
        z = z/z0
        alpha1 = 0
        alpha3 = 1
        dif_vec = x - alpha3*z
        g3 = evalg(dif_vec)

        while g3>=g1:
            alpha3 = alpha3/2
            dif_vec = x - alpha3*z
            g3 = evalg(dif_vec)
            
        if alpha3<tol:
            print("no likely improvement")
            ier = 0
            return [x,g1,ier]
        
        alpha2 = alpha3/2
        dif_vec = x - alpha2*z
        g2 = evalg(dif_vec)

        h1 = (g2 - g1)/alpha2
        h2 = (g3-g2)/(alpha3-alpha2)
        h3 = (h2-h1)/alpha3

        alpha0 = 0.5*(alpha2 - h1/h3)
        dif_vec = x - alpha0*z
        g0 = evalg(dif_vec)

        if g0<=g3:
            alpha = alpha0
            gval = g0

        else:
            alpha = alpha3
            gval =g3

        x = x - alpha*z

        if abs(gval - g1)<tol:
            ier = 0
            return [x,gval,ier]

    print('max iterations exceeded')    
    ier = 1        
    return [x,g1,ier]

def driver2():

    print ('For Newton:')
    t = time.time()
    for j in range(50):
        [xstar1,ier1,its1] =  Newton(x,tol, Nmax)
    elapsed1 = time.time()-t
    print(xstar1)
    print('Newton: the error message reads:',ier1) 
    print('Newton: took this many seconds:',elapsed1/50)
    print('Netwon: number of iterations is:',its1)


    print('For Steepest Descent:')
    t = time.time()
    for j in range(50):
        [xstar,gval,ier] =  SteepestDescent(x,tol,Nmax)
    elapsed1 = time.time()-t
    print("the steepest descent code found the solution ",xstar)
    print('Newton: took this many seconds:',elapsed1/50)
    print("g evaluated at this point is ", gval)
    print("ier is ", ier)


    tolnew = 5*10**-2
    xnew = SteepestDescent(x,tolnew,Nmax)[0]
    print(xnew)


    print('For Steepest Descent Part 2:')
    t = time.time()
    for j in range(50):
        [xstar,gval,ier] =  SteepestDescent(xnew,tol,Nmax)
    elapsed1 = time.time()-t
    print("the steepest descent code found the solution ",xstar)
    print('Newton: took this many seconds:',elapsed1/50)
    print("g evaluated at this point is ", gval)
    print("ier is ", ier)

driver2()

