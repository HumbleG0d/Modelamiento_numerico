import numpy as np


def cholesky_decomposition(A, b):
    # Verificamos si la matriz es defninida positiva
    minors_nonzero = all(
        np.linalg.det(A[: i + 1, : i + 1]) > 0 for i in range(A.shape[0])
    )
    # Verficamos si la matriz es simétrica
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
        print("Matriz triangular inferior L:")
        print(L)

        # Verificar que L * L^T sea igual a A
        print("Verificación: L * L^T == A")
        print(np.allclose(L @ L.T, A))

        """
      Tenemos el sistema Ax=b.
      Mediante cholesKy tenemos que A=LL^T
      Entonces nuestro nuevo sistema es LL^Tx = b
      Sea L^Tx = y. Entonces tenemos los siguientes sistemas
      Ly=b ; L^Tx = y
      Primero resolvamos Ly = b , luego L^Tx = y.
    """
        print("\nCalculando la matriz ´y´ por sustitucion progresiva.")
        y = np.zeros(n)
        for i in range(n):
            y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]

        print("Matriz y:")
        print(y)

        # Finlmente resolvemos L^Tx = y. Para ello aplicamos sustitucion regresiva
        print("\nCalculando la matriz ´x´ por sustitucion regresiva.")
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - np.dot(L.T[i, i + 1 :], x[i + 1 :])) / L.T[i, i]

        print("Matriz solucion:")
        print(x)

    else:
        print("no se puede resolver por Cholosky")


# Matriz de ejemplo
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
b = np.array(
    [
        [56],
        [25],
        [62],
        [41],
        [10],
        [47],
    ],
    dtype=float,
)

# Obtener la descomposición de Cholesky
cholesky_decomposition(A, b)
