import numpy as np


def jacobi_richardson(A, b, x0, max_iterations, tolerance):
    n = len(b)
    x = x0.copy()

    for iteration in range(max_iterations):
        x_new = np.zeros(n)

        for i in range(n):
            sigma = 0
            for j in range(n):
                if i != j:
                    sigma += A[i, j] * x[j]
            x_new[i] = (b[i] - sigma) / A[i, i]
        print(f"Iteracion {iteration + 1} : ", x_new)
        if np.linalg.norm(x_new - x) < tolerance:
            print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
            return x_new

        x = x_new

    print("Convergencia no alcanzada después de", max_iterations, "iteraciones.")
    return x


# Ejemplo de uso:
A = np.array([[10.0, 1.0], [2.0, 10.0]])
b = np.array([11.0, 12.0])
x0 = np.zeros_like(b, dtype=float)  # Aproximación inicial como float
max_iterations = 100  # Número máximo de iteraciones
tolerance = 1e-6  # Tolerancia de convergencia

solution = jacobi_richardson(A, b, x0, max_iterations, tolerance)
print("Solución aproximada:", solution)
