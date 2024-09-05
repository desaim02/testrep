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
    plt.legend()  # Show legend to differentiate between plots
    plt.show()

driver()

# The difference is that (x-2)^9 with polynomials is much less precise than the original (x-2)^9 function. 
# The discrepency is caused because of the repeated subtractions which is losing accuracy. 
# In the original function, there is only subtraction being done. 
# Whereas the expanded function has 5 subtractions being done and with each subtraction, precision is lost. 
# The "most correct" plot would be (x-2)^9. 


# Problem 3 

import numpy as np

def f(x):
    return (1+x+x**3)*(np.cos(x))

def P(x):
    return 1+x-((x**2)/2)
   #return  1 + x*(((1+x+x**3)*(-np.sin(x))) + (np.cos(x)*(1+3*x**2))) + ((.5*x**2)*((x**3)+6*x+x-1)*np.cos(x))

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

y2og = [original(x2,S) for S in Delta]
y2manip = [manip(x2,S) for S in Delta]

plt.plot (Delta,y1og, label= "original function")
plt.plot (Delta,y1manip, label= 'manipulated function')
plt.xscale('log')
plt.xlabel('S')
plt.ylabel('y')
plt.title(x = {x1})
plt.legend()  
plt.show()

plt.plot (Delta,y2og, label= "original function")
plt.plot (Delta,y2manip, label= 'manipulated function')
plt.xscale('log')
plt.xlabel('S')
plt.ylabel('y')
plt.title(x = {x2})
plt.legend()  
plt.show()


# C

def algo(x,S):
    return -S*np.sin(x) + .5*S**2*np.cos(E) #PUT E IN

y1algo = [algo(x1,S) for S in Delta]
y2algo = [algo(x2,S) for S in Delta]

plt.plot (Delta,y1og, label= "original function")
plt.plot (Delta,y1manip, label= 'manipulated function')
plt.plot (Delta,y1algo,label="algorithm")
plt.xscale('log')
plt.xlabel('S')
plt.ylabel('y')
plt.title(x = {x1})
plt.legend()  
plt.show()

plt.plot (Delta,y2og, label= "original function")
plt.plot (Delta,y2manip, label= 'manipulated function')
plt.plot (Delta,y2algo,label="algorithm")
plt.xscale('log')
plt.xlabel('S')
plt.ylabel('y')
plt.title(x = {x2})
plt.legend()  
plt.show()



# Part b

# Function to compute direct difference and modified expression
def cos_diff_direct(x, delta):
    return np.cos(x + delta) - np.cos(x)

def cos_diff_manip(x, delta):
    return -2 * np.sin((2*x + delta) / 2) * np.sin(delta / 2)

# Define x values
x_values = [np.pi, 1e6]

# Define delta values
delta_values = np.logspace(-16, 0, num=16)

# Plotting
plt.figure(figsize=(10, 5))

for x in x_values:
    direct_diff = [cos_diff_direct(x, delta) for delta in delta_values]
    manip_diff = [cos_diff_manip(x, delta) for delta in delta_values]

    plt.loglog(delta_values, np.abs(np.array(direct_diff) - np.array(manip_diff)), label=f'x = {x}')
    
plt.xlabel('Delta')
plt.ylabel('Absolute Difference')
plt.title('Difference between Direct Subtraction and Manipulated Expression')
plt.legend()

plt.savefig("HW1.5.b.png")
plt.clf

# Part c

# Function using Taylor expansion
def cos_diff_taylor(x, delta):
    return -delta * np.sin(x) + (delta ** 2 / 2) * (-np.cos(x))

# Compare Taylor expansion with other methods
plt.figure(figsize=(10, 5))

for x in x_values:
    taylor_diff = [cos_diff_taylor(x, delta) for delta in delta_values]
    manip_diff = [cos_diff_manip(x, delta) for delta in delta_values]

    plt.loglog(delta_values, np.abs(np.array(taylor_diff) - np.array(manip_diff)), label=f'x = {x}')

plt.xlabel('Delta')
plt.ylabel('Absolute Difference')
plt.title('Difference between Taylor Expansion and Manipulated Expression')
plt.legend()

plt.savefig("HW1.5.c.png")

