import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad

def driver():
    
    f = lambda s: 1/(1+(s**2))          
    a = -5
    b = 5
    
    # exact integral
    I_ex = 2*math.atan(5)
    
#    N =100
#    ntest = np.arrange(0,N,step=2)
    
#    errorT = np.zeros(len(ntest))
#    errorS = np.zeros(len(ntest))
    
#    for j in range(0,len(ntest)):
#        n = ntest[j]

# for simpson's n must be even.        
# n+1 = number of pts.
    n = 220
    I_trap = CompTrap(a,b,n,f)
    print('I_trap= ', I_trap)
    
    err = abs(I_ex-I_trap)   
    
    print('absolute error = ', err)    

    print('relative error = ', err/I_trap)
    
    I_simp = CompSimp(a,b,n,f)

    print('I_simp= ', I_simp)
    
    err = abs(I_ex-I_simp)   
    
    print('absolute error = ', err)    
        
def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    
    I_trap = h*f(xnode[0])*1/2
    
    for j in range(1,n):
         I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    
    return I_trap     

def CompSimp(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_simp = f(xnode[0])

    nhalf = n/2
    for j in range(1,int(nhalf)+1):
         # even part 
         I_simp = I_simp+2*f(xnode[2*j])
         # odd part
         I_simp = I_simp +4*f(xnode[2*j-1])
    I_simp= I_simp + f(xnode[n])
    
    I_simp = h/3*I_simp
    
    return I_simp   

# 2

def driver2():
    #t = 1/x
    f = lambda t: t * np.cos(1/t)
    a = 10**-16 #close to 0 but not exact
    b = 1
    n = 4
    I_ex = 0.01811762198060567 #from integal calculator
    
    I_simp = CompSimp(a,b,n,f)
    print('I_simp= ', I_simp)
    err = abs(I_ex-I_simp)   
    print('absolute error = ', err) 

driver2()