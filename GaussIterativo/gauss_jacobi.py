import numpy as np


def gauss_jacobi(A, b, x0, max_iterations, tolerance):
    n = len(b)
    x = x0.copy()

    for iteration in range(max_iterations):
        x_new = np.zeros(n)

        for i in range(n):
            sigma = 0
            for j in range(n):
                if i != j:
                    sigma += A[i, j] * x[j]
            x_new[i] = (1 / A[i, i]) * (b[i] - sigma)

        error = np.linalg.norm(x_new - x)
        print(
            f"Iteración {iteration + 1}: Solución aproximada = {x_new}, Error = {error}"
        )

        if error < tolerance:
            print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
            return x_new

        x = x_new

    print("Convergencia no alcanzada después de", max_iterations, "iteraciones.")
    return x


def jacobi_matrix(A):
    n = A.shape[0]
    D = np.diag(np.diag(A))  # Extraer la matriz diagonal de A
    L = np.tril(A, k=-1)  # Extraer la matriz triangular inferior
    U = np.triu(A, k=1)  # Extraer la matriz triangular superior

    # Calcular la matriz de Jacobi
    J = -np.linalg.inv(D).dot(L + U)

    return J


# Ejemplo de uso:
A = np.array(
    [
        [4.0, -1.0, 0.0, -1.0, 0.0, 0.0],
        [-1.0, 4.0, -1.0, 0.0, -1.0, 0.0],
        [0.0, -1.0, 4.0, 0.0, 0.0, -1.0],
        [-1.0, 0.0, 0.0, 4.0, -1.0, 0.0],
        [0.0, -1.0, 0.0, -1.0, 4.0, -1.0],
        [0.0, 0.0, -1.0, 0.0, -1.0, 4.0],
    ]
)

b = np.array([56.0, 25.0, 62.0, 41.0, 10.0, 47.0])

x0 = np.zeros_like(b, dtype=float)  # Aproximación inicial como float
max_iterations = 100  # Número máximo de iteraciones
tolerance = 1e-6  # Tolerancia de convergencia

solution = gauss_jacobi(A, b, x0, max_iterations, tolerance)
J = jacobi_matrix(A)
print("Matriz de Jacobi:")
print(J)
print("Solución aproximada:", solution)
