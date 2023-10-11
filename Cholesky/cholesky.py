import numpy as np


def cholesky_decomposition(A):
    n = A.shape[0]
    L = np.zeros_like(A, dtype=float)

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                L[i, j] = np.sqrt(A[i, i] - np.sum(L[i, :] ** 2))
            else:
                L[i, j] = (A[i, j] - np.sum(L[i, :] * L[j, :])) / L[j, j]

    return L


# Matriz de ejemplo
A = np.array([[3, 1, 21], [1, 9, 41], [21, 41, 401]], dtype=float)

# Obtener la descomposición de Cholesky
L = cholesky_decomposition(A)

# Mostrar la matriz triangular inferior L
print("Matriz triangular inferior L:")
print(L)

# Verificar que L * L^T sea igual a A
print("Verificación: L * L^T == A")
print(np.allclose(L @ L.T, A))
