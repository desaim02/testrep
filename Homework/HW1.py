import numpy as np 
import matplotlib.pyplot as plt

#Problem 1 

def poly (x):
    return (x**9) - (18*(x**8)) + (144*(x**7)) - (672*(x**6)) + (2016*(x**5)) - (4032*(x**4)) + (5376*(x**3)) - (4608*(x**2)) + (2304*x) -512
 
def driver():
    xvals = np.arange(1.920,2.080,.001)
    y1 = [poly(x) for x in xvals]
    y2 = [(x-2)**9 for x in xvals]

    plt.plot(xvals, y1, label='with polynomials')
    plt.plot(xvals, y2, label='original', linestyle='--')  

    # Label the axes
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('(x-2)^9')
    plt.legend()  

driver()

# Problem 3 

import numpy as np

def f(x):
    return (1+x+x**3)*(np.cos(x))

def P(x):
    return 1+x-((x**2)/2)

print("P(.5) is", P(.5))
print("f(.5) is", f(.5))

print(" The actual error is", abs(f(.5) - P(.5)))

#Problem 4

a = 1
b = -56
c = 1

discriminant = b**2 -4*a*c
discrimround = round(np.sqrt(discriminant),3)
discrimnr = np.sqrt(discriminant)

grround = (-b + discrimround)/2
brround = (-b - discrimround)/2
print(grround,brround)

grnr= (-b + discrimnr)/2
brnr = (-b - discrimnr)/2
print(grnr,brnr)


abserrgr = abs(grnr-grround)
abserrbr = abs(brnr-brround)
print(abserrgr,abserrbr)


relerrgr = abserrgr/grnr
relerrbr = abserrbr/brnr
print(relerrgr,relerrbr)


op1 =  -grround-b
op2 = c/grround
print(op1,op2)

ab1 = abs(brnr-op1)
ab2 = abs(brnr-op2)
print(ab1,ab2)

rel1 = ab1/brnr
rel2 = ab2/brnr
print(rel1,rel2)


# Problem 5

# A


# B
import numpy as np
import matplotlib.pyplot as plt

def original(x,S):
    return np.cos(x+S) - np.cos(x)

def manip(x,S):
    return -2*np.sin((2*x+S)/2)*np.sin(S/2)

x1 = np.pi
x2 = 10**6

DeltaExp = np.arange(-16,1,1)
Delta = 10**DeltaExp.astype(float)

y1og = [original(x1,S) for S in Delta]
y1manip = [manip(x1,S) for S in Delta]
diff1 = np.array(y1og)-np.array(y1manip)

y2og = [original(x2,S) for S in Delta]
y2manip = [manip(x2,S) for S in Delta]
diff2 = np.array(y2og) -np.array(y2manip)

plt.plot(Delta,diff1)
plt.xscale('log')
plt.title(f"Difference between original and manipulated function for x = pi")
plt.show()

plt.plot(Delta,diff2)
plt.xscale('log')
plt.title(f"Difference between original and manipulated function for x = {x2}")
plt.show()

# C

def algo(x,S):
    return -S*np.sin(x)

y1algo = [algo(x1,S) for S in Delta]
y2algo = [algo(x2,S) for S in Delta]

diff1pt2 = np.array(y1og) - np.array(y1algo)
diff2pt2 = np.array(y2og) - np.array(y2algo)

# plt.plot(Delta,diff1pt2)
# plt.xscale('log')
# plt.title(f"Difference between original and algorithm for x = pi")
# plt.show()

# plt.plot(Delta,diff2pt2)
# plt.xscale('log')
# plt.title(f"Difference between original and algorithm for x = {x2}")
# plt.show()

print(Delta)
print(diff1)
print(diff1pt2)

print(Delta)
print(diff2)
print(diff2pt2)

