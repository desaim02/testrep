import numpy as np


def newton(f,fp,p0,tol,Nmax):
  
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1)
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def newton_matrix(f, g, jacobian, x0, y0, tol=1e-13, nmax=100):
    x, y = x0, y0
    for n in range(nmax):
        F = np.array([f(x, y), g(x, y)])
        J = jacobian(x, y)

        # Check if the Jacobian is singular (to avoid division by zero)
        if np.linalg.det(J) == 0:
            raise ValueError("Jacobian is singular, cannot invert.")

        delta = np.linalg.solve(J, -F)

        x += delta[0]
        y += delta[1]


        if np.linalg.norm(delta) < tol:
            print(f"Converged after {n + 1} iterations.")
            return x, y

    print("Did not converge within the maximum number of iterations.")
    return x, y


# Question 1 


def iteration(f, g, M, x0, y0, tol=1e-13, nmax=100):
    x, y = x0, y0
    n = 0
    diff = tol + 1  
    iterate_list = []

    while diff > tol and n < nmax:
        f_val = f(x, y)
        g_val = g(x, y)
        F = np.array([f_val, g_val])
        xy_new = np.array([x, y]) - M @ F
        iterate_list.append(xy_new)
        diff = np.linalg.norm(xy_new - np.array([x, y]))
        x, y = xy_new[0], xy_new[1]
        n += 1
    
    if diff <= tol:
        print(f"Converged after {n} iterations.")
    else:
        print("Did not converge within the maximum number of iterations.")
    
    print(iterate_list)
    return x, y

def driver1a():
    x0, y0 = 1.0, 1.0
    f = lambda x, y: 3 * x**2 - y**2
    g = lambda x, y: 3 * x * y**2 - x**3 - 1
    M = np.array([[1/6, 1/18],
                [0, 1/6]])

    solution_x, solution_y = iteration(f, g, M, x0, y0)

    print(f"Solution: x = {solution_x}, y = {solution_y}")

driver1a()



#Question 1 c
f = lambda x, y: 3 * x**2 - y**2
g = lambda x, y: 3 * x * y**2 - x**3 - 1
jacobian = lambda x, y: np.array([[6 * x, -2 * y],
                                    [3 * y**2 - 3 * x**2, 6 * x * y]])
x0, y0 = 1.0, 1.0

solution_x, solution_y = newton_matrix(f, g, jacobian, x0, y0)

print(f"Solution: x = {solution_x}, y = {solution_y}")


# Question 3c

def iteration(f, fpx, fpy, fpz, d, x0, y0, z0, tol=1e-13, Nmax=100):
    n = 0
    vec0 = np.array([x0, y0, z0])
    results = [vec0]  

    iter_x_values = np.array([[n, x0]])  

    while n < Nmax:
        x1 = x0 - d(x0, y0, z0) * fpx(x0, y0, z0)
        y1 = y0 - d(x0, y0, z0) * fpy(x0, y0, z0)
        z1 = z0 - d(x0, y0, z0) * fpz(x0, y0, z0)

        vec1 = np.array([x1, y1, z1])
        results.append(vec1)  

        iter_x_values = np.vstack([iter_x_values, [n + 1, x1]])

        delta = vec1 - vec0

        if np.linalg.norm(delta) < tol:
            print(f"Converged after {n + 1} iterations.")
            return vec1, np.array(results), iter_x_values

        vec0 = vec1
        x0, y0, z0 = x1, y1, z1

        n += 1

    print("Did not converge.")
    return vec1, np.array(results), iter_x_values


# Define the functions
f = lambda x, y, z: x**2 + 4*y**2 + 4*z**2 - 16
fpx = lambda x, y, z: 2*x
fpy = lambda x, y, z: 8*y
fpz = lambda x, y, z: 8*z

# Initial values
x0 = y0 = z0 = 1

# Define d with dynamic values
d = lambda x, y, z: f(x, y, z) / (fpx(x, y, z)**2 + fpy(x, y, z)**2 + fpz(x, y, z)**2)

# Call the iteration function
final_result, results, iter_x_values = iteration(f, fpx, fpy, fpz, d, x0, y0, z0)

print(final_result, results)





