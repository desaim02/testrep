import numpy as np

# Problem 1 

def driver1a():
    f = lambda x,y,z: 6*x + 2*y + 2*z + 2
    g = lambda x,y,z: 2*x + (2/3)*y+(1/3)*z-1
    h = lambda x,y,z: x + 2*y-z

    x,y,z = 2.6,-3.8,-5

    # if this is 0 (or close to 0), then solution is verified
    print(f(x,y,z),g(x,y,z),h(x,y,z))

#driver1a()

# Problem 3

# a

# Power Method
def power(A,tol,nmax):
    # choose a random x_0 based on matrix A dimensions 
    x_0 = np.random.rand(A.shape[1])
    it_ct = 0
    eval_actual = np.max(np.linalg.eigvals(A))
    x_1_norm = 100
    while abs(x_1_norm - eval_actual) > tol and it_ct < nmax:
            it_ct += 1
            x_1 = np.dot(A,x_0)
            x_1_norm = np.linalg.norm(x_1)
            # normalized vector becomes new inital
            x_0 = x_1/x_1_norm
    return x_1_norm,x_0,it_ct

# Create Hilbert Matrices
def Hilbert_create(n):
    mat = np.zeros((n, n))
    for i in range(1, n+1):
        for j in range(1, n+1):
            mat[i-1, j-1] = 1 / (i + j - 1)
    return mat

def driver3a():
    tol = 10**-15
    nmax = 100

    for n in range(4,21,4):
        print('for n ', n )
        mat = Hilbert_create(n)
        eval_actual = np.max(np.linalg.eigvals(mat))
        eval_power,evec_power,it_ct = power(mat,tol,nmax)
        #print("actual eval is", evec_power)
        #print("actual evec is", eval_actual)
        #print("power estimate is", eval_power)
        print("error is", eval_actual-eval_power)
        print("iteration count is", it_ct)

#driver3a()

#3b

#Implement Inverse Power Method
def inverse_power(A, tol, nmax):
    # Choose a random x_0 based on matrix A dimensions
    x_0 = np.random.rand(A.shape[1])
    it_ct = 0
    eval_old = 0
    x_1_norm = 100  # Start with a dummy large value for the eigenvalue approximation

    while abs(x_1_norm - eval_old) > tol and it_ct < nmax:
        it_ct += 1
        # Solve A*y = x_0 for y instead of explicitly computing A^-1
        y = np.linalg.solve(A, x_0)
        x_1_norm = np.linalg.norm(y)  # Dominant eigenvalue of A^-1 (smallest eigenvalue of A)
        eval_old = x_1_norm
        x_0 = y / x_1_norm  # Normalize the vector

    smallest_eigenvalue = 1 / x_1_norm  # Inverse gives the smallest eigenvalue
    return smallest_eigenvalue, x_0, it_ct

def driver_3b():
    tol = 10**-8  # Convergence tolerance
    nmax = 100  # Maximum number of iterations
    n = 16  # Matrix size

    mat = Hilbert_create(n)
    eval_power,evec_power,it_ct = inverse_power(mat, tol, nmax)

    # Compare with actual eigenvalues
    eval_actual = np.min(np.linalg.eigvals(mat))

    # Display results
    print("actual eval is", eval_actual)
    #print("actual evec is", eval_actual)
    print("power estimate is", eval_power)
    print("error is", eval_actual-eval_power)
    print("iteration count is", it_ct)


# Run the driver function
#driver_3b()

def driver_3b_error_analysis():
    tol = 1e-8  
    nmax = 100  
    n = 16  

    mat = Hilbert_create(n)
    eval_power, evec_power, it_ct = inverse_power(mat, tol, nmax)
    eval_actual = np.min(np.linalg.eigvals(mat))

    # Estimate the perturbation induced by numerical errors
    perturbed_mat = mat.copy()
    perturbed_mat[0, 0] += 1e-12  # Simulate a tiny perturbation
    perturbed_eval = np.min(np.linalg.eigvals(perturbed_mat))
    perturbation_norm = np.linalg.norm(perturbed_mat - mat, ord=2)

    actual_error = abs(eval_actual - eval_power)

    print("Actual smallest eigenvalue:", eval_actual)
    print("Power method estimate:", eval_power)
    print("Error (|λ_actual - λ_computed|):", actual_error)
    print("Simulated perturbation eigenvalue shift:", abs(eval_actual - perturbed_eval))
    print("Perturbation norm (∥E∥):", perturbation_norm)
    print("Bauer-Fike bound:", perturbation_norm)
    print("Iteration count:", it_ct)

# Run the function
driver_3b_error_analysis()
