import numpy as np


def cholesky_decomposition(A):
    minors_nonzero = all(
        np.linalg.det(A[: i + 1, : i + 1]) > 0 for i in range(A.shape[0])
    )
    transpuesta = np.array_equal(A, A.T)

    if minors_nonzero and transpuesta:
        n = A.shape[0]
        L = np.zeros_like(A, dtype=float)

        for i in range(n):
            for j in range(i + 1):
                if i == j:
                    L[i, j] = np.sqrt(A[i, i] - np.sum(L[i, :] ** 2))
                else:
                    L[i, j] = (A[i, j] - np.sum(L[i, :] * L[j, :])) / L[j, j]

        # Mostrar la matriz triangular inferior L

        # Verificar que L * L^T sea igual a A
        return L

    else:
        print("no se puede resolver por Cholosky")


def LDLT(A):
    """
    Function to perform the LDL^T matrix decomposition
    """
    C = cholesky_decomposition(A)
    n = len(A)
    D = np.identity(n, dtype=np.float64)
    L = C
    for i in range(0, n):
        D[i, i] = C[i, i] ** 2
        L[:, i] = C[:, i] * D[i, i] ** (-1 / 2)
    return L, D


A = np.array(
    [
        [4, -1, 0, -1, 0, 0],
        [-1, 4, -1, 0, -1, 0],
        [0, -1, 4, 0, 0, -1],
        [-1, 0, 0, 4, -1, 0],
        [0, -1, 0, -1, 4, -1],
        [0, 0, -1, 0, -1, 4],
    ],
    dtype=float,
)


L, D = LDLT(A)
print("Matriz L:")
print(L)
print("Matriz D:")
print(D)
