import numpy as np 
import matplotlib.pyplot as plt
import random
import math

#Question 2 

#b 

A = .5*np.array([[1, 1],
              [1 + 10**-10, 1 - 10**-10]])

Ainv = np.array([[1-10**10, 10**10],
                       [1+10**10, - 10**10]])

norm_A_2 = np.linalg.norm(A, ord=2)
norm_Ainv_2 = np.linalg.norm(Ainv, ord=2)

cond_num = norm_A_2*norm_Ainv_2

print(cond_num)

#c 

b = np.array([[1],
              [1]])

print(np.linalg.norm(b, ord=2))


#Question 3

#c

def y(x):
    return math.exp(x)

def f(x):
    return y(x)-1

x1 = 9.999999995000000 * 10**-10
realx = 10**-9

errrel = (realx - f(x1))/realx
print(errrel)

# xvals = np.linspace(0,100,100)
# range(len(xvals))

# yvals = [y(x) for x in xvals]
# fvals = [f(x) for x in xvals]

# plt.plot(xvals,yvals)
# plt.plot(xvals,fvals)
# plt.show()

# diff = [yvals[i]-fvals[i] for i in list(range(len(xvals)))]

# plt.plot(xvals,diff)
# plt.show()

#d


def taylore(x,n):
    result = 0
    for i in range(1, n + 1):
        result += (x ** i) / math.factorial(i)
    return result


x1 = 9.999999995000000 * 10**-10
realx = 10**-9

n = 1
errrel = 1  

while errrel > 10**-16:
    taylor_result = taylore(x1, n)  
    errrel = abs(realx - taylor_result) / abs(realx)  
    print("n = ", n, "errel = ", errrel)  
    n += 1
    
print(n-1) # have to subtract one because the n is incremented at the end of the while loop

print(taylore(x1,2))


#Question 4

#a
t = np.arange(0,np.pi,np.pi/30)
y = [np.cos(tentry) for tentry in t]

N = 3 # insert N 
kvec = np.arange(1,N,1)
sumvec = [t[k]*y[k] for k in kvec]
sumtot = np.sum(sumvec)
print("the sum is: ", sumtot)

#b
R = 1.2
Sr = .1
f = 15
p = 0

def x(theta):
    return R*(1+Sr*np.sin(f*theta+p))*np.cos(theta)

def y(theta):
    return R*(1+Sr*np.sin(f*theta+p))*np.sin(theta)


thetavals = np.linspace(0,2*np.pi,100)

xvals = [x(theta) for theta in thetavals]
yvals = [y(theta) for theta in thetavals]

plt.plot(thetavals,xvals)
plt.plot(thetavals,yvals)
plt.xlim(0, 6)  
plt.ylim(-10,10) 
plt.show()


ivec = np.arange(0,10,1)
Sr = .05
Rvec = [i for i in ivec]
fvec = [2*i for i in ivec]


xresult = [[Rvec[i]*(1+Sr*np.sin(fvec[i]*theta+p))*np.cos(theta) for i in ivec] for theta in thetavals]
yresult = [[Rvec[i]*(1+Sr*np.sin(fvec[i]*theta+p))*np.sin(theta) for i in ivec] for theta in thetavals]

plt.plot(thetavals,xresult)
plt.plot(thetavals,yresult)
plt.xlim(0,6)  
plt.ylim(-10, 10)  
plt.show()


print(((10**20)*2+(10**10)*2+1)/(10**20))
print(2+(10**-10))