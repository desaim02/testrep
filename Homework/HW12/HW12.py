import numpy as np

# Problem 1 

#a
f = lambda x,y,z: 6*x + 2*y + 2*z + 2
g = lambda x,y,z: 2*x + (2/3)*y+(1/3)*z-1
h = lambda x,y,z: x + 2*y-z

x,y,z = 2.6,-3.8,-5

# if this is 0 (or close to 0), then solution is verified
print(f(x,y,z),g(x,y,z),h(x,y,z))

