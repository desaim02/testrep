import numpy as np 
import matplotlib.pyplot as plt
import random

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