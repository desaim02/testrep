import matplotlib.pyplot as plt
import numpy as np
# PreLab

'''
Write a subroutine that constructs and evaluates a line that goes through the points (x0, f (x0)) and
(x1, f (x1)) at a point α.
'''

def line_eval(x0,x1,y0,y1):
    ydiff = y1-y0
    xdiff = x1-x0
    m = ydiff/xdiff
    line = lambda x: m*(x-x0)+y0
    #print('the line is', line)
    
    # xvals = np.linspace(x0,x1,100)
    # yvals = f(xvals)
    # plt.plot(xvals,yvals, label = 'actual')
    # plt.plot(x0,y0,'o',color='red')
    # plt.plot(x1,y1,'o',color='red')
    # plt.plot(xvals,line(xvals), label = 'approx')


    # plt.legend()
    # plt.show()

    return line


def  eval_lin_spline(xeval,Neval,a,b,f,Nint):

    '''create the intervals for piecewise approximations'''
    xint = np.linspace(a,b,Nint+1)
   
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval) 

    for jint in range(Nint):
        ind = xint[jint]

        if xeval
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        '''let n denote the length of ind'''
        
        '''temporarily store your info for creating a line in the interval of 
         interest'''
         a1= xint[jint]
         fa1 = f(a1)
         b1 = xint[jint+1]
         fb1 = f(b1)
        
        for kk in range(n):
           #yeval(ind(kk)) = line_eval(xeval(ind(kk)))
           '''use your line evaluator to evaluate the lines at each of the points 
           in the interval'''
           '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
           the points (a1,fa1) and (b1,fb1)'''
           
           
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()               


def drivertest():
    f = lambda x: x**3-2
    x0 = 1.88
    x1 = 0
    y0 = f(x0)
    y1= f(x1)

    
    line = line_eval(x0,x1,y0,y1)
    print(line(x0))
    print(y0)

drivertest()



def exercise32()
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    # Nint = 
    # xeval = 
    # Neval = 

    #eval_lin_spline(xeval,Neval,a,b,f,Nint)
