import numpy as np


def sor(A, b, x0, omega, max_iterations, tolerance):
    n = len(b)
    x = x0.copy()

    for iteration in range(max_iterations):
        for i in range(n):
            suma = np.dot(A[i, :], x) - A[i, i] * x[i]
            x[i] = (1 - omega) * x[i] + (omega / A[i, i]) * (b[i] - suma)

        error = np.linalg.norm(np.dot(A, x) - b)
        print(f"Iteración {iteration + 1}: Solución aproximada = {x}, Error = {error}")

        if error < tolerance:
            print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
            return x

    print("Convergencia no alcanzada después de", max_iterations, "iteraciones.")
    return x


# Ejemplo de uso:
A = np.array([[4.0, 3.0, 0.0], [3.0, 4.0, -1.0], [0.0, -1.0, 4.0]])
b = np.array([24.0, 30.0, -24.0])
x0 = np.zeros_like(b, dtype=float)  # Aproximación inicial como float
omega = 1.25  # Parámetro de relajación
max_iterations = 100  # Número máximo de iteraciones
tolerance = 1e-6  # Tolerancia de convergencia

solution = sor(A, b, x0, omega, max_iterations, tolerance)
print("Solución aproximada:", solution)
