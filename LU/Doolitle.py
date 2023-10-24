import numpy as np


def lu_doolitle(A, b):
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
            L[i][i] = 1
            print("l", i + 1, i + 1, "=", L[i][i])
            if i == 0:
                U[0][0] = A[0][0]
                for j in range(1, n):
                    U[0][j] = A[0][j]
                    L[j][0] = A[j][0] / U[0][0]
                    print("l", j + 1, 1, "=", L[j][0])
                    print("u", 1, j + 1, "=", U[0][j])
            else:
                for j in range(i, n):
                    temp = 0
                    for k in range(0, i):
                        temp = temp + L[i][k] * U[k][j]
                    U[i][j] = A[i][j] - temp
                    print("u", i + 1, j + 1, "=", U[i][j])

                for j in range(i + 1, n):
                    temp = 0
                    for k in range(0, i):
                        temp = temp + L[j][k] * U[k][i]
                    L[j][i] = (A[j][i] - temp) / U[i][i]
                    print("l", j + 1, i + 1, "=", L[j][i])

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


A = np.array(
    [
        [625, 125, 25, 5, 1],
        [256, 64, 16, 4, 1],
        [81, 27, 9, 3, 1],
        [16, 8, 4, 2, 1],
        [1, 1, 1, 1, 1],
    ],
    dtype=float,
)

b = np.array([[225], [100], [36], [9], [1]], dtype=float)

lu_doolitle(A, b)
