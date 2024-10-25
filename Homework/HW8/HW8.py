import numpy as np
import math
from numpy.linalg import inv 
from numpy.linalg import norm
import matplotlib.pyplot as plt

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

def eval_hermite(xeval,xint,yint,ypint,N):


    ''' Evaluate all Lagrange polynomials'''

    lj = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    ''' Construct the l_j'(x_j)'''
    lpj = np.zeros(N+1)
#    lpj2 = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
#              lpj2[count] = lpj2[count]*(xint[count] - xint[jj])
              lpj[count] = lpj[count]+ 1./(xint[count] - xint[jj])
              

    yeval = 0.
    
    for jj in range(N+1):
       Qj = (1.-2.*(xeval-xint[jj])*lpj[jj])*lj[jj]**2
       Rj = (xeval-xint[jj])*lj[jj]**2
#       if (jj == 0):
#         print(Qj)
         
#         print(Rj)
#         print(Qj)
#         print(xeval)
 #        return
       yeval = yeval + yint[jj]*Qj+ypint[jj]*Rj
       
    return(yeval)
    
def create_clamped_spline(yint,yintp,xint,N): #added yintp which is the derivatives at all the points

#    create the right hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)  
    for i in range(1,N): 
       hi = xint[i]-xint[i-1]
       hip = xint[i+1] - xint[i]
       b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
       h[i-1] = hi
       h[i] = hip

#   since most of the rows are the same, keep this for loop and change the first and last row of b vector
    #h0 = xint[1]-xint[0]
    #hn_1 = xint[N]-xint[N-1]
    b[0] = -yintp[0] + ((yint[1]-yint[0])/h[0])
    b[N] = -yintp[N] + ((yint[N]-yint[N-1])/h[N-1])


#  create matrix so you can solve for the M values
# This is made by filling one row at a time 
    A = np.zeros((N+1,N+1))
    A[0][0] = h[0]/3
    A[0][1] = h[0]/6
    for j in range(1,N):
       A[j][j-1] = h[j-1]/6
       A[j][j] = (h[j]+h[j-1])/3 
       A[j][j+1] = h[j]/6
    A[N][N] = h[N-1]/3
    A[N][N-1] = h[N-1]/6
    #print(A[0],A[N])

    Ainv = inv(A)
    
    M  = Ainv.dot(b)

#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = yint[j]/h[j]-h[j]*M[j]/6
       D[j] = yint[j+1]/h[j]-h[j]*M[j+1]/6
    return(M,C,D)
       
def create_natural_spline(yint,xint,N):

#    create the right hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)  
    for i in range(1,N):
       hi = xint[i]-xint[i-1]
       hip = xint[i+1] - xint[i]
       b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
       h[i-1] = hi
       h[i] = hip

#  create matrix so you can solve for the M values
# This is made by filling one row at a time 
    A = np.zeros((N+1,N+1))
    A[0][0] = 1.0
    for j in range(1,N):
       A[j][j-1] = h[j-1]/6
       A[j][j] = (h[j]+h[j-1])/3 
       A[j][j+1] = h[j]/6
    A[N][N] = 1

    Ainv = inv(A)
    
    M  = Ainv.dot(b)

#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = yint[j]/h[j]-h[j]*M[j]/6
       D[j] = yint[j+1]/h[j]-h[j]*M[j+1]/6
    return(M,C,D)
       
def eval_local_spline(xeval,xi,xip,Mi,Mip,C,D):
# Evaluates the local spline as defined in class
# xip = x_{i+1}; xi = x_i
# Mip = M_{i+1}; Mi = M_i

    hi = xip-xi
    yeval = (Mi*(xip-xeval)**3 +(xeval-xi)**3*Mip)/(6*hi) \
            + C*(xip-xeval) + D*(xeval-xi)
    return yeval 
        
def  eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    
    yeval = np.zeros(Neval+1)
    
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        
#   find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]

# evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
#        print('yloc = ', yloc)
#   copy into yeval
        yeval[ind] = yloc

    return(yeval)    

# DRIVERS 

#def driver_es_LH():
    f = lambda x: 1./(1.+x**2)
    fp = lambda x: -2*x/(1.+x**2)**2

    n = [5,10,15,20]

    ''' interval'''
    a = -5
    b = 5
   
    for N in n:
        ''' create equispaced interpolation nodes'''
        xint = np.linspace(a,b,N+1)
        
        ''' create interpolation data'''
        yint = np.zeros(N+1)
        ypint = np.zeros(N+1)
        for jj in range(N+1):
            yint[jj] = f(xint[jj])
            ypint[jj] = fp(xint[jj])
        
        ''' create points for evaluating the Lagrange interpolating polynomial'''
        Neval = 1000
        xeval = np.linspace(a,b,Neval+1)
        yevalL = np.zeros(Neval+1)
        yevalH = np.zeros(Neval+1)
        for kk in range(Neval+1):
            yevalL[kk] = eval_lagrange(xeval[kk],xint,yint,N)
            yevalH[kk] = eval_hermite(xeval[kk],xint,yint,ypint,N)

        ''' create vector with exact values'''
        fex = np.zeros(Neval+1)
        for kk in range(Neval+1):
            fex[kk] = f(xeval[kk])
        
        
        plt.figure()
        plt.plot(xeval,fex,'ro-',label="Actual")
        plt.plot(xeval,yevalL,'bs--',label='Lagrange') 
        plt.plot(xeval,yevalH,'c.--',label='Hermite')
        plt.title(f"N={N}")
        plt.legend()
        plt.semilogy()
        plt.show()
            
        errL = abs(yevalL-fex)
        errH = abs(yevalH-fex)
        plt.figure()
        plt.semilogy(xeval,errL,'bs--',label='Lagrange')
        plt.semilogy(xeval,errH,'c.--',label='Hermite')
        plt.legend()
        plt.title(f"N={N}")
        plt.show() 


#EQUISPACED DRIVERS
def driver_es_LH():
    f = lambda x: 1./(1.+x**2)
    fp = lambda x: -2*x/(1.+x**2)**2
    n = [5, 10, 15, 20]
    a, b = -5, 5

    fig, axes = plt.subplots(4, 2, figsize=(16, 8))  # Create a 2x4 grid for 8 plots
    fig.suptitle("Lagrange and Hermite for Equispaced Nodes")  # Overall title

    for idx, N in enumerate(n):
        xint = np.linspace(a, b, N + 1)
        yint = f(xint)
        ypint = fp(xint)
        
        Neval = 1000
        xeval = np.linspace(a, b, Neval + 1)
        yevalL = np.array([eval_lagrange(x, xint, yint, N) for x in xeval])
        yevalH = np.array([eval_hermite(x, xint, yint, ypint, N) for x in xeval])
        
        fex = f(xeval)
        
        # Plot actual function vs. Lagrange and Hermite interpolations
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label="Actual")
        ax1.plot(xeval, yevalL, 'bs--', label='Lagrange')
        ax1.plot(xeval, yevalH, 'c.--', label='Hermite')
        ax1.set_title(f"N={N}")
        ax1.legend()
        
        # Plot error in semilog scale
        errL = abs(yevalL - fex)
        errH = abs(yevalH - fex)
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, errL, 'bs--', label='Lagrange')
        ax2.semilogy(xeval, errH, 'c.--', label='Hermite')
        ax2.set_title(f"Error, N={N}")
        ax2.legend()

    plt.tight_layout()
    plt.show()


def driver_es_cubic_natural():
    f = lambda x: 1./(1.+x**2)
    a = -5
    b = 5

    n = [5, 10, 15, 20]
    
    fig, axes = plt.subplots(4, 2, figsize=(16, 8))  # Create a 2x4 grid for 8 plots
    fig.suptitle("Natural Cubic for Equispaced Nodes")  # Overall title

    for idx, Nint in enumerate(n):
        xint = np.linspace(a, b, Nint + 1)
        yint = f(xint)

        # Create points to evaluate at
        Neval = 100
        xeval = np.linspace(xint[0], xint[Nint], Neval + 1)

        (M, C, D) = create_natural_spline(yint, xint, Nint)

        yeval = eval_cubic_spline(xeval, Neval, xint, Nint, M, C, D)

        # Evaluate f at the evaluation points
        fex = f(xeval)
        nerr = norm(fex - yeval)

        # Plot exact function vs natural spline
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label='Exact function')
        ax1.plot(xeval, yeval, 'bs--', label='Natural spline') 
        ax1.legend()
        ax1.set_title(f"N={Nint}")

        # Plot absolute error in semilog scale
        err = abs(yeval - fex)
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, err, 'ro--', label='Absolute error')
        ax2.legend()
        ax2.set_title(f"Error, N={Nint}")

    plt.tight_layout()  # Adjust layout
    plt.show()


def driver_es_cubic_clamped():
    f = lambda x: 1/(1+x**2)
    fp = lambda x: -2*x/(1+x**2)**2
    a = -5
    b = 5

    n = [5, 10, 15, 20]

    fig, axes = plt.subplots(4, 2, figsize=(16, 8))  # Create a 2x4 grid for 8 plots
    fig.suptitle("Clamped Cubic for Equispaced Nodes")  # Overall title

    for idx, Nint in enumerate(n):
        xint = np.linspace(a, b, Nint + 1)
        yint = f(xint)
        yintp = fp(xint)

        # Create points to evaluate at
        Neval = 100
        xeval = np.linspace(xint[0], xint[Nint], Neval + 1)

        (M, C, D) = create_clamped_spline(yint, yintp, xint, Nint)

        yeval = eval_cubic_spline(xeval, Neval, xint, Nint, M, C, D)

        # Evaluate f at the evaluation points
        fex = f(xeval)
        nerr = norm(fex - yeval)

        # Plot exact function vs cubic spline
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label='Exact function')
        ax1.plot(xeval, yeval, 'bs--', label='Clamped Cubic spline')
        ax1.legend()
        ax1.set_title(f"N={Nint}")

        # Plot absolute error in semilog scale
        err = abs(yeval - fex)
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, err, 'ro--', label='Absolute error')
        ax2.legend()
        ax2.set_title(f"Error, N={Nint}")

    plt.tight_layout()  # Adjust layout
    plt.show()

#CHEBYSHEV DRIVERS  

def driver_ch_LH():
    f = lambda x: 1. / (1. + x ** 2)
    fp = lambda x: -2 * x / (1. + x ** 2) ** 2

    n = [5, 10, 15, 20]

    # Create a 2x4 grid for 8 plots
    fig, axes = plt.subplots(4, 2, figsize=(16, 8))
    fig.suptitle("Lagrange and Hermite for Chebyshev Nodes")  # Overall title

    a = -5
    b = 5

    for idx, N in enumerate(n):
        # Create Chebyshev interpolation nodes
        xint = [5 * np.cos(((2 * k + 1) * np.pi) / (2 * (N + 1))) for k in range(N + 1)]
        print(xint)

        # Create interpolation data
        yint = np.zeros(N + 1)
        ypint = np.zeros(N + 1)
        for jj in range(N + 1):
            yint[jj] = f(xint[jj])
            ypint[jj] = fp(xint[jj])

        # Create points for evaluating the Lagrange interpolating polynomial
        Neval = 1000
        xeval = np.linspace(a, b, Neval + 1)
        yevalL = np.zeros(Neval + 1)
        yevalH = np.zeros(Neval + 1)
        for kk in range(Neval + 1):
            yevalL[kk] = eval_lagrange(xeval[kk], xint, yint, N)
            yevalH[kk] = eval_hermite(xeval[kk], xint, yint, ypint, N)

        # Create vector with exact values
        fex = np.zeros(Neval + 1)
        for kk in range(Neval + 1):
            fex[kk] = f(xeval[kk])

        # Plot exact function vs Lagrange and Hermite
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label="Actual")
        ax1.plot(xeval, yevalL, 'bs--', label='Lagrange')
        ax1.plot(xeval, yevalH, 'c.--', label='Hermite')
        ax1.set_title(f"N={N}")
        ax1.legend()
        
        # Calculate errors
        errL = abs(yevalL - fex)
        errH = abs(yevalH - fex)

        # Plot absolute error in semilog scale
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, errL, 'bs--', label='Lagrange')
        ax2.semilogy(xeval, errH, 'c.--', label='Hermite')
        ax2.set_title(f"Error, N={N}")
        ax2.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout for overall title
    plt.show()


def driver_ch_cubic_natural():
    f = lambda x: 1. / (1. + x ** 2)
    a = -5
    b = 5

    n = [5, 10, 15, 20]

    # Create a 2x4 grid for plots
    fig, axes = plt.subplots(4, 2, figsize=(16, 8))
    fig.suptitle("Natural Cubic Spline Interpolation")  # Overall title

    for idx, Nint in enumerate(n):
        xint = np.array([5 * np.cos(((2 * k + 1) * np.pi) / (2 * (Nint + 1))) for k in range(Nint + 1)])
        xint = np.sort(xint)
        print(xint)
        yint = f(xint)

        # Create points for evaluation
        Neval = 100
        xeval = np.linspace(xint[0], xint[Nint], Neval + 1)

        (M, C, D) = create_natural_spline(yint, xint, Nint)

        print('M =', M)

        yeval = eval_cubic_spline(xeval, Neval, xint, Nint, M, C, D)

        # Evaluate f at the evaluation points
        fex = f(xeval)

        nerr = norm(fex - yeval)
        print('nerr =', nerr)

        # Plot exact function vs natural spline
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label='Exact Function')
        ax1.plot(xeval, yeval, 'bs--', label='Natural Spline')
        ax1.set_title(f"N={Nint}")
        ax1.legend()

        # Calculate error
        err = abs(yeval - fex)

        # Plot absolute error in semilog scale
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, err, 'ro--', label='Absolute Error')
        ax2.set_title(f"Error, N={Nint}")
        ax2.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout for overall title
    plt.show()


def driver_ch_cubic_clamped():
    f = lambda x: 1 / (1 + x ** 2)
    fp = lambda x: -2 * x / (1 + x ** 2) ** 2
    a = -5
    b = 5

    n = [5, 10, 15, 20]

    # Create a 2x4 grid for plots
    fig, axes = plt.subplots(4, 2, figsize=(16, 8))
    fig.suptitle("Clamped Cubic Spline Interpolation")  # Overall title

    for idx, Nint in enumerate(n):
        xint = [5 * np.cos(((2 * k + 1) * np.pi) / (2 * (Nint + 1))) for k in range(Nint + 1)]
        xint = np.sort(xint)
        yint = [f(x) for x in xint]
        yintp = [fp(x) for x in xint]

        # Create points for evaluation
        Neval = 100
        xeval = np.linspace(xint[0], xint[Nint], Neval + 1)

        (M, C, D) = create_clamped_spline(yint, yintp, xint, Nint)

        print('M =', M)

        yeval = eval_cubic_spline(xeval, Neval, xint, Nint, M, C, D)

        # Evaluate f at the evaluation points
        fex = f(xeval)

        nerr = norm(fex - yeval)
        print('nerr =', nerr)

        # Plot exact function vs clamped cubic spline
        ax1 = axes.flat[2 * idx]
        ax1.plot(xeval, fex, 'ro-', label='Exact Function')
        ax1.plot(xeval, yeval, 'bs--', label='Cubic Spline')
        ax1.set_title(f"N={Nint}")
        ax1.legend()

        # Calculate error
        err = abs(yeval - fex)

        # Plot absolute error in semilog scale
        ax2 = axes.flat[2 * idx + 1]
        ax2.semilogy(xeval, err, 'ro--', label='Absolute Error')
        ax2.set_title(f"Error, N={Nint}")
        ax2.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout for overall title
    plt.show()

# DRIVER CALLS
# driver_es_LH()
# driver_es_cubic_natural()
# driver_es_cubic_clamped()
# driver_ch_LH()
#driver_ch_cubic_natural()
#driver_ch_cubic_clamped()