import numpy as np 
import matplotlib.pyplot as plt

x = [1,2,3]
3*x

#When you multiple a vector by 3, it duplicates the elements in a vectro 3 times. 

y=np.array([1,2,3])
3*y

#Multiplying an array by 3 multiplies all the elements within the array by 3. 

print('this is 3y',3*y)

X = np.linspace(0,2 * np.pi,100)
Ya = np.sin(X)
Yb = np.cos(X)

plt.plot(X,Ya)
plt.plot(X,Yb)
#plt.show()

print(X.size)

#The size of X is 100. np.linspace(x,y,z) breaks up the range between x and y with z points. 

plt.xlabel('x')
plt.ylabel('y')
plt.show()

## Exercises 

#Exercise 1 

x = np.linspace(0,10,10)
y = np.arange(0,10,1)

print(len(x))
print(len(y))

#Exercise 2

print(x[:3])

#Exercise 3

print('the first three entries of x are', x[:3])

#Exercise 4

w = 10**(-np.linspace(1,10,10))
print(w)
xtoo = np.arange(0,len(w),1)

plt.semilogy(xtoo,w)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#Exercise 5 

s = 3*w
plt.semilogy(xtoo,w)
plt.semilogy(xtoo,s)
plt.xlabel('x')
plt.ylabel('y')

#Exercise 6 

plt.savefig('Exercise5')
plt.show()

#Practical Code Design 
#Sample Code
import numpy as np
import numpy.linalg as la
import math

def dotProduct(x,y,n):
    # Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp

def driver():
    n = 100
    x = np.linspace(0,np.pi,n)
    # this is a function handle. You can use it to define
    # functions instead of using a subroutine like you
    # have to in a true low level language.
    
    f = lambda x: x**2 + 4*x + 2*np.exp(x)
    g = lambda x: 6*x**3 + 2*np.sin(x)
    y = f(x)
    print(y)
    w = g(x)
    print(w)
    # evaluate the dot product of y and w
    dp = dotProduct(y,w,n)
    # print the output
    print('the dot product is : ', dp)
    return

driver()

#4.2.1
def dotProduct(x,y,n):
    # Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp

vector1 = np.array([1, 0])
vector2 = np.array([0, 1])
n = len(vector1)
result = dotProduct(vector1, vector2, n)
print(result)

#4.2.2
def dotProduct(x,y,n):
    # Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp

random_matrix1 = np.random.rand(3, 2)
random_matrix2 = np.random.rand(2, 3)
print(random_matrix1)
print(random_matrix2)

def matrixmult(x,y):
    if x.shape[1] != y.shape[0]:
        raise ValueError("Matrices are not aligned for multiplication")
    n = x.shape[1]
    result = np.zeros((x.shape[0], y.shape[1]))
    for i in range(x.shape[0]):          
        for j in range(y.shape[1]):     
            result[i, j] = dotProduct(x[i, :], y[:, j],n)
    return result

def driver():   
   mm = matrixmult(random_matrix1,random_matrix2)
   print(mm)

driver()

#4.2.3
random_vector1 = np.random.rand(3)
random_vector2 = np.random.rand(3)
np.dot(random_vector1,random_vector2)
np.dot(random_matrix1,random_matrix2)



#Python Background Knowledge
x = np.zeros((2,3))
print(x)
print ('len(x):', len(x))
print ('x.size:', x.size)
print ('x.shape:', x.shape)
x[:,1]







