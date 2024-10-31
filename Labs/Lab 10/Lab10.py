import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy.linalg as la
import math

def eval_legendre(n,x):
    leg_return = np.zeros(n+1)
    leg_return[0] = 1
    if n>0:
        leg_return[1] = x
        for i in range(1,n):
            leg1 = leg_return[i]
            legm1 = leg_return[i-1]
            #f = lambda x: (1/(i+1))* (2*i+1)*x*leg1(x) - (n*legm1(x))
            leg_return[i+1] = ((2*i+1)*x*leg1- (i*legm1)) / (i+1)
    print(leg_return)
    return leg_return

#eval_legendre(3,1)

def coeff(f,phi,w,a,b):
    return quad(phi*f*w,a,b) / quad(phi**2*w,a,b)

def driver():
#  function you want to approximate
    f = lambda x: math.exp(x)

# Interval of interest    
    a = -1
    b = 1
# weight function    
    w = lambda x: 1.

# order of approximation
    n = 2

#  Number of points you want to sample in [a,b]
    N = 100
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])
      
    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])
        
    plt.figure()    
    plt.plot(xeval,fex,'ro-', label= 'f(x)')
    plt.plot(xeval,pval,'bs--',label= 'Expansion') 
    plt.legend()
    plt.show()    
    
    err = abs(pval-fex)
    plt.semilogy(xeval,err_l,'ro--',label='error')
    plt.legend()
    plt.show()
    
      
def eval_legendre_expansion(f,a,b,w,n,x): 

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab 
  p = eval_legendre(n,x) 
  # initialize the sum to 0 
  pval = 0.0    
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_legendre(n,x)[j]  #already evaluated at x within function 
      # make a function handle for evaluating phi_j^2(x)*w(x)
      phi_j_sq = lambda x: eval_legendre(n,x)[j]**2 * w(x)
      # use the quad function from scipy to evaluate normalizations
      norm_fac,err = quad(phi_j_sq,a,b)
      # make a function handle for phi_j(x)*f(x)*w(x)/norm_fac
      func_j = lambda x: (phi_j(x)*f(x)*w(x))/norm_fac
      # use the quad function from scipy to evaluate coeffs
      aj,err = quad(func_j,a,b)
      # accumulate into pval
      pval = pval+aj*p[j] 
       
  return pval


def eval_legendre_expansion2(f,a,b,w,n,x): 

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab 
  p = eval_legendre(n,x) 
  # initialize the sum to 0 
  pval = 0.0   
  a = np.zeros(n+1) 
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_legendre(n,x)[j]  #already evaluated at x within function 
      # make a function handle for evaluating phi_j^2(x)*w(x)
      phi_j_sq = lambda x: eval_legendre(n,x)[j]**2 * w(x)
      # use the quad function from scipy to evaluate normalizations
      norm_fac,err = quad(phi_j_sq,a,b)
      # make a function handle for phi_j(x)*f(x)*w(x)/norm_fac
      func_j = lambda x: (phi_j(x)*f(x)*w(x))/norm_fac
      # use the quad function from scipy to evaluate coeffs
      a[j],err = quad(func_j,a,b)
      # accumulate into pval

  for j in range(0,n+1):
        pval = pval+a[j]*p[j] 
       
  return pval
    

driver()