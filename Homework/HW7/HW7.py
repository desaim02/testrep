import numpy as np 
import matplotlib.pyplot as plt
import numpy.linalg as la
from numpy.linalg import inv
from numpy.linalg import norm

def  eval_monomial(xeval,coef,N,Neval):

    yeval = coef[0]*np.ones(Neval+1)
    
#    print('yeval = ', yeval)
    
    for j in range(1,N+1):
      for i in range(Neval+1):
#        print('yeval[i] = ', yeval[i])
#        print('a[j] = ', a[j])
#        print('i = ', i)
#        print('xeval[i] = ', xeval[i])
        yeval[i] = yeval[i] + coef[j]*xeval[i]**j

    return yeval
  
  
def Vandermonde(xint,N):

    V = np.zeros((N+1,N+1))
    
    ''' fill the first column'''
    for j in range(N+1):
       V[j][0] = 1.0

    for i in range(1,N+1):
        for j in range(N+1):
           V[j][i] = xint[j]**i

    return V     

# 1b 

# Code for seeing multiple plots on one 
def drivermult(N,ax): 

    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    
    ''' Create interpolation nodes'''
    h = 2/(N-1)

    xint = np.zeros(N+1)
    for i in range (1, N+2):
        xint[i-1] =  -1 + (i-1)*h

   # xint = np.linspace(a,b,N+1)
#    print('xint =',xint)
    '''Create interpolation data'''
    yint = f(xint)
#    print('yint =',yint)
    
    ''' Create the Vandermonde matrix'''
    V = Vandermonde(xint,N)
#    print('V = ',V)

    ''' Invert the Vandermonde matrix'''    
    Vinv = inv(V)
#    print('Vinv = ' , Vinv)
    
    ''' Apply inverse to rhs'''
    ''' to create the coefficients'''
    coef = Vinv @ yint
    
#    print('coef = ', coef)

# No validate the code
    Neval = 100    
    xeval = np.linspace(a,b,Neval+1)
    yeval = eval_monomial(xeval,coef,N,Neval)


#plot function
    # plt.plot(xint,f(xint),'o',label=f'Nodes N={N}')
    # plt.plot(xeval,yeval,label=f'Interpolation N={N}') #coefficient functon 
    # plt.plot(xeval,f(xeval), '--', label=f'Exact function N={N}')

    #plt.title("X = {N}")
    #plt.show()

# exact function
    # yex = f(xeval)
    
    # err =  norm(yex-yeval) 
    # print('err = ', err)

    ax.plot(xint, yint, 'o', label=f'Nodes N={N}')
    ax.plot(xeval, yeval, label=f'Interpolation N={N}')
    ax.plot(xeval, f(xeval), '--', label=f'Exact function N={N}')
    ax.set_title(f"N = {N}")
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    #ax.legend(loc="upper center")
    return

def drivermult2():
    N=  range(11,20)
    fig, axes = plt.subplots(1, len(N), figsize=(15, 5))  # 1 row, len(Ns) columns
    # Loop over each N and its corresponding subplot axis
    for i, N in enumerate(N):
        drivermult(N, axes[i])  # Pass each axis (axes[i]) to the driver function

    # Adjust layout so plots don't overlap
    plt.tight_layout()
    plt.legend()
    plt.show()

drivermult2()


def driver1(N): 

    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    
    ''' Create interpolation nodes'''
    #confused on the x_i and h business (doesn't it just end up being on the interval -1 and 1 )
    #h = 2/N-1
    #x_i = lambda i,h: -1 + (i-1)*h
    xint = np.linspace(a,b,N+1)
#    print('xint =',xint)
    '''Create interpolation data'''
    yint = f(xint)
#    print('yint =',yint)
    
    ''' Create the Vandermonde matrix'''
    V = Vandermonde(xint,N)
#    print('V = ',V)

    ''' Invert the Vandermonde matrix'''    
    Vinv = inv(V)
#    print('Vinv = ' , Vinv)
    
    ''' Apply inverse to rhs'''
    ''' to create the coefficients'''
    coef = Vinv @ yint
    
#    print('coef = ', coef)

# No validate the code
    Neval = 1000    
    xeval = np.linspace(a,b,Neval+1)
    yeval = eval_monomial(xeval,coef,N,Neval)

#plot function
    # plt.plot(xint,f(xint),'o',label=f'Nodes N={N}')
    # plt.plot(xeval,yeval,label=f'Interpolation N={N}') #coefficient functon 
    # plt.plot(xeval,f(xeval), '--', label=f'Exact function N={N}')
    # plt.legend()
    # plt.show()

# exact function
    # yex = f(xeval)
    
    # err =  norm(yex-yeval) 
    # print('err = ', err)

    print(coef)

    return

#driver1(20)

#2 

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)

''' create divided difference matrix'''
def dividedDiffTable(x, y, n):
 
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /
                                     (x[j] - x[i + j]));
    return y;
    
def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)
    
    ptmp[0] = 1.
    for j in range(N):
      ptmp[j+1] = ptmp[j]*(xval-xint[j])
     
    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N+1):
       yeval = yeval + y[0][j]*ptmp[j]  

    return yeval

def lagrange_bary(x_eval, x_interp, y_interp, degree, func):

    # temp_val = np.ones(degree+1)
    
    # for index in range(degree+1):
    #    for j in range(degree+1):
    #        if (j != index):
    #           temp_val[index] = temp_val[index]*(x_eval - x_interp[j])/(x_interp[index]-x_interp[j])

    # y_result = 0
    
    # for j in range(degree+1):
    #    y_result = y_result + y_interp[j]*temp_val[j]

    product_term = 1

    for x_int in x_interp:
        if x_int != x_eval:
            product_term *= (x_eval - x_int)

    total_sum = 0

    for j in range(degree+1):

        denom = 1

        for i in range(degree+1):
            if x_interp[j] != x_interp[i]:
                denom *= (x_interp[j] - x_interp[i])
            
        weight = 1 / denom

        total_sum += (weight / (x_eval - x_interp[j])) * func(x_interp[j])
    
    y_result = product_term * total_sum

    return y_result



def driver2(N):
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
   
    xint = np.linspace(a,b,N+1)
    yint = f(xint)
    

    Neval = 1000
    xeval = np.linspace(a,b,Neval+1)
    yeval_l= np.zeros(Neval+1)
    yeval_lb= np.zeros(Neval+1)
    #yeval_dd = np.zeros(Neval+1)
  
    '''Initialize and populate the first columns of the 
     divided difference matrix. We will pass the x vector'''
   # y = np.zeros( (N+1, N+1) )
     
    # for j in range(N+1):
    #    y[j][0]  = yint[j]

    # y = dividedDiffTable(xint, y, N+1)
    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
       #yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
       yeval_lb[kk] = lagrange_bary(xeval[kk],xint,yint,N,f)
       #yeval_dd[kk] = evalDDpoly(xeval[kk],xint,y,N)
          

    ''' create vector with exact values'''
    fex = f(xeval)
       

    plt.figure()    
    plt.plot(xint,yint,'o',label=f'Nodes N={N}',markersize=10, markeredgewidth=2, markeredgecolor='black', markerfacecolor='red')
    plt.plot(xeval,fex,'ro-', label = "exact function")
    #plt.plot(xeval,yeval_l,'bs--', label = "Lagrange Approx") 
    plt.plot(xeval,yeval_lb,'bs--', label = "Lagrange Barry Approx") 
    #plt.plot(xeval,yeval_dd,'c.--', label = "Divided Diff Approx")
    plt.legend()

    # plt.figure() 
    # err_l = abs(yeval_l-fex)
    # err_dd = abs(yeval_dd-fex)
    # plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    # plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
    # plt.legend()
    plt.show()
    return 

driver2(17)
driver2(18)
driver2(19)

#3 
# use monomial basis (part a)
def driver3(N):
    xint = np.zeros(N+1)
    for i in range(1,N+1):
       xint[i-1] = np.cos((((2*i)-1)*np.pi)/(2*N))
       
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1

    '''Create interpolation data'''
    yint = f(xint)
    #    print('yint =',yint)

    ''' Create the Vandermonde matrix'''
    V = Vandermonde(xint,N)
    #    print('V = ',V)

    ''' Invert the Vandermonde matrix'''    
    Vinv = inv(V)
    #    print('Vinv = ' , Vinv)

    ''' Apply inverse to rhs'''
    ''' to create the coefficients'''
    coef = Vinv @ yint

    #    print('coef = ', coef)
    # No validate the code
    Neval = 1000  
    xeval = np.linspace(a,b,Neval+1)
    yeval = eval_monomial(xeval,coef,N,Neval)

    # exact function
    yex = f(xeval)

    plt.plot(xint,yint,'o',label=f'Nodes N={N}')
    plt.plot(xeval,yeval,label=f'Interpolation N={N}') #coefficient functon 
    plt.plot(xeval,yex, '--', label=f'Exact function N={N}')
    plt.legend()
    plt.title("X = {N}")
    plt.show()

    # err =  norm(yex-yeval) 
    # print('err = ', err)

    return

driver3(17)
driver3(18)


# Trying to find why even vs odd N's produce such different results. 
# def xfunc(N):
#     f = lambda x: 1/(1+(10*x)**2)
#     xarr =np.zeros(N+1)
#     for i in range(1,N+1,2):
#         xarr[N-1] = np.cos((((2*i)-1)*np.pi)/(2*N))
#     plt.plot(xarr,f(xarr))
#     plt.show()
    
# xfunc(18)

