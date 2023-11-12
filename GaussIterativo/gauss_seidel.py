import numpy as np


def gauss_seidel(A, b, x0, max_iterations, tolerance):
    n = len(b)
    x = x0.copy()

    for iteration in range(max_iterations):
        for i in range(n):
            suma = np.dot(A[i, :], x) - A[i, i] * x[i]
            x[i] = (1 / A[i, i]) * (b[i] - suma)

        error = np.linalg.norm(np.dot(A, x) - b)
        print(f"Iteración {iteration + 1}: Solución aproximada = {x}, Error = {error}")

        if error < tolerance:
            print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
            return x

    print("Convergencia no alcanzada después de", max_iterations, "iteraciones.")
    return x


def gauss_seidel_matrix(A):
    n = A.shape[0]
    D = np.diag(np.diag(A))  # Extraer la matriz diagonal de A
    L = np.tril(A, k=-1)  # Extraer la matriz triangular inferior
    U = np.triu(A, k=1)  # Extraer la matriz triangular superior

    # Calcular la matriz de Gauss-Seidel
    GS = -np.dot(np.linalg.inv(D - L), U)

    return GS


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


solution = gauss_seidel(A, b, x0, max_iterations, tolerance)
print("Solución aproximada:", solution)
GS = gauss_seidel_matrix(A)
print("Matriz de Gauss-Seidel:")
print(GS)
