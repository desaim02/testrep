import numpy as np
import pandas as pd
import time
import rpy2.robjects as ro

'''
I think everyhting is set up rn, but results aren't good because covariance matrix is highly ill-conditioned for both Boston and Iris.
Will explore to see if I can find better dataset.
'''

# Load dataset
ro.r('library(MASS)')
ro.r('data <- readRDS("/Users/mitalidesai/Desktop/APPM 4600/ProjectCode/iris.rds")')
r_iris = ro.r('data')
iris_df = pd.DataFrame.from_dict(r_iris)
X = iris_df

# Direct PCA with eigenvalue decomposition
def direct_pca(X):
    # Step 1: Center the data (zero mean)
    X_centered = X - X.mean(axis=0)

    # Step 2: Compute the covariance matrix
    cov_matrix = np.cov(X_centered, rowvar=False)
    print ('condition number for covariance matrix:', np.linalg.cond(cov_matrix))

    # Step 3: Measure time for eigendecomposition
    start_time = time.time()
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    end_time = time.time()

    # Time taken for direct eigendecomposition
    time_taken = end_time - start_time
    print(f"Time taken for direct PCA (eigendecomposition): {time_taken:.6f} seconds")

    # Sort eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues_sorted = eigenvalues[sorted_indices]

    return eigenvalues_sorted, time_taken


# Arnoldi iteration function
def arnoldi_iteration(A, v0, m):
    n = A.shape[0]
    Q = np.zeros((n, m))
    H = np.zeros((m, m))
    Q[:, 0] = v0 / np.linalg.norm(v0)

    for k in range(m - 1):
        v = np.dot(A, Q[:, k])
        for j in range(k + 1):
            H[j, k] = np.dot(Q[:, j], v)
            v = v - H[j, k] * Q[:, j]

        H[k, k] = np.linalg.norm(v)
        if H[k, k] > 1e-10:
            Q[:, k + 1] = v / H[k, k]
        else:
            break
    return Q, H

def arnoldi_eigenvalues(A, v0, m):
    Q, H = arnoldi_iteration(A, v0, m)
    eigenvalues = np.linalg.eigvals(H)
    return eigenvalues

def pca_with_arnoldi(A, m, v0=None):
    A_centered = A - np.mean(A, axis=0)
    cov_matrix = np.cov(A_centered, rowvar=False)

    if v0 is None:
        v0 = np.random.rand(cov_matrix.shape[0])

    eigenvalues = arnoldi_eigenvalues(cov_matrix, v0, m)

    # Sort eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues_sorted = eigenvalues[sorted_indices]

    return eigenvalues_sorted


# Set the number of principal components (eigenvalues)
m = 5

# 1. Run PCA with direct eigenvalue computation
direct_eigenvalues, direct_time = direct_pca(X)

# 2. Run PCA with Arnoldi iteration
start_time = time.time()
arnoldi_eigenvalues_sorted = pca_with_arnoldi(X, m)
arnoldi_time = time.time() - start_time

# Output the results
print("\nDirect PCA Eigenvalues:")
print(direct_eigenvalues[:5])
print(f"Time taken for direct PCA: {direct_time:.6f} seconds\n")

print("\nPCA with Arnoldi Eigenvalues:")
print(arnoldi_eigenvalues_sorted[:5])
print(f"Time taken for Arnoldi PCA: {arnoldi_time:.6f} seconds")

