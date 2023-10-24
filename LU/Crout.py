import numpy as np


def lu_crout(A, b):
    # Verficamos que la determinante de cada matriz esquina sea diferente de cero
    minors_nonzero = all(
        np.linalg.det(A[: i + 1, : i + 1]) != 0 for i in range(A.shape[0])
    )

    if not minors_nonzero:
        print("La matriz no se puede descomponer por LU")
    else:
        m, n = A.shape
        L = np.zeros((n, n))
        U = np.zeros((n, n))
        s1 = 0
        s2 = 0
        for i in range(n):
            L[i][0] = A[i][0]
            U[i][i] = 1
            print("U", i + 1, i + 1, "=", U[i][i])
            print("L", i + 1, 1, "=", L[i][0])

        for j in range(1, n):
            U[0][j] = A[0][j] / L[0][0]
            print("U", 1, j + 1, "=", U[0][j])
        for k in range(1, n):
            for i in range(k, n):
                for r in range(k):
                    s1 += L[i][r] * U[r][k]
                L[i][k] = A[i][k] - s1
                print("L", i + 1, k + 1, "=", L[i][k])
                s1 = 0
            for j in range(k + 1, n):
                for r in range(k):
                    s2 += L[k][r] * U[r][j]
                U[k][j] = (A[k][j] - s2) / L[k][k]
                print("U", k + 1, j + 1, "=", U[k][j])
                s2 = 0
        print("\nMatriz L:")
        print(L)
        print("\nMatriz U:")
        print(U)

        """
      Inicialmenete teniamos el sistema Ax = b
      Ahora tenemos el sistema LUx = b.
      Sea Ux = y , asi tendremos Ly = b
      Primero resolvamos Ly = b.
      Luego resolvemos Ux = y.
    """
        # Resolviendo Ly = b.Para ello aplicamos sustitucion progresiva
        print("\nCalculando la matriz ´y´ por sustitucion progresiva.")
        y = np.zeros(n)
        for i in range(n):
            y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]

        print("Matriz y:")
        print(y)

    # Finlmente resolvemos Ux = y. Para ello aplicamos sustitucion regresiva
    print("\nCalculando la matriz ´x´ por sustitucion regresiva.")
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i + 1 :], x[i + 1 :])) / U[i, i]

    print("Matriz solucion:")
    print(x)


A = np.array([[10, -3, 6], [1, -8, -2], [-2, 4, 9]], dtype=float)

b = np.array([[25], [-9], [-50]], dtype=float)

lu_crout(A, b)
