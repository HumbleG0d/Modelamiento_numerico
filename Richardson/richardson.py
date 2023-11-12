import numpy as np


def richardson(A, b, x0, alpha, max_iterations, tolerance):
    x = x0
    for iteration in range(max_iterations):
        Ax = np.dot(A, x)
        residual = Ax - b
        x -= alpha * residual
        print(f"Iteracion {iteration + 1} : ", x)
        # Verificar la condición de convergencia
        if np.linalg.norm(residual - x0) < tolerance:
            print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
            break
    # print(f"Convergencia alcanzada después de {iteration + 1} iteraciones.")
    return x


# Ejemplo de uso:
A = np.array([[0.5, -0.2, 0.5], [0.1, 0.6, 0.4], [-0.3, 0.1, 0.0]])
b = np.array([-1.0, 6.5, 0.7])
x0 = np.zeros_like(b)  # Aproximación inicial
alpha = 0.1  # Tamaño de paso
max_iterations = 100  # Número máximo de iteraciones
tolerance = 1e-6  # Tolerancia de convergencia

solution = richardson(A, b, x0, alpha, max_iterations, tolerance)
print("Solución aproximada:", solution)
