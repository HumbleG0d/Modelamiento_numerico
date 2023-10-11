import numpy as np


def lu_doolittle(A):
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros((n, n))

    for i in range(n):
        # Calcular U
        for j in range(i, n):
            U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])

        # Calcular L
        for j in range(i + 1, n):
            L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]

    return L, U


A = np.array([[1, 2, 20], [1, -2, 0], [1, 1, 1]], dtype=float)


L, U = lu_doolittle(A)


print("Matriz L (triangular inferior):")
print(L)

print("\nMatriz U (triangular superior):")
print(U)


print("\nVerificación: A == L * U")
print(np.allclose(A, L @ U))  # Verficamos
