import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad


def driver1():
    xvals = np.linspace(0,5,100)
    Maclaurin = lambda x: x - x**3/6 +x**5/math.factorial(5) 
    acpade = lambda x: (x - (7/60)*x**3) / (1+(1/20)*x**2)
    bpade = lambda x: x / (1+(1/6)*x**2+(7/360)*x**4)
    realfunc = lambda x: np.sin(x)


    maclaurin_err = Maclaurin(xvals) - realfunc(xvals)
    acpade_err = acpade(xvals) - realfunc(xvals)
    bpade_err = bpade(xvals) - realfunc(xvals)


    #plt.plot(xvals,realfunc(xvals), label = "actual function")
    # Plot Error Curves
    plt.plot(xvals,maclaurin_err, label = "Maclaurin")
    plt.plot(xvals,acpade_err, label = "Parts A and C Pade")
    plt.plot(xvals,bpade_err, label = "Part B Pade")
    plt.legend()
    plt.title('Error of Pade and Maclaurin for f(x) = sin(x)')
    plt.show()

    # Plot Actual Curves
    plt.plot(xvals, realfunc(xvals), label = 'Actual')
    plt.plot(xvals,Maclaurin(xvals), label = "Maclaurin")
    plt.plot(xvals,acpade(xvals), label = "Parts A and C Pade")
    plt.plot(xvals,bpade(xvals), label = "Part B Pade")
    plt.legend()
    plt.title('Pade and Maclaurin for f(x) = sin(x)')
    plt.show()



def quad_l_Trap(f,a,b):
    h = (b-a)/2
    xvals = [a,b]
    fvec = np.array([f(x) for x in xvals])
    w = np.array([h/2,h/2])
    return sum(w*fvec)

def quad_l_Simp(f,a,b):
    h = (b-a)/2
    xvals = [a,(a+b)/2,b]
    fvec = np.array([f(x) for x in xvals])
    w = np.array([h/3,4*h/3,h/3])
    return sum(w*fvec)

def driver3a(type):
    a = -5
    b = 5
    n = 10
    f = lambda s: 1/(1+(s**2))          
    xvals = np.linspace(a,b,n)
    exact = 2*math.atan(5)
    to_sum = []
    for i in range(n-1):
        if type == "Trap":
            to_sum.append(quad_l_Trap(f,xvals[i],xvals[i+1]))
        if type == "Simp":
            to_sum.append(quad_l_Simp(f,xvals[i],xvals[i+1]))
    print(f'for type {type} the sum is {sum(to_sum)}')
    print(f' the absoutle error is {sum(to_sum)-exact}')

    return sum(to_sum)





def driver3c():
    a = -5
    b = 5
    f = lambda s: 1/(1+(s**2))         
    result, error,infodict = quad(f,a,b,full_output=True)
    neval = infodict['neval']
    print (f' for the default tolerance, the result is {result}, the error is {error}, neval: {neval}')
    result2, error2,infodict2 = quad(f,a,b,epsrel = 10**-4,full_output=True)
    neval2 = infodict2['neval']
    print (f' for the modified tolerance, the result is {result2}, the error is {error2}, neval: {neval2}')

#driver1()
#driver3a('Trap')
#driver3a('Simp')
#driver3c()



#PROF CODE - used for 3c
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

driver()
    


