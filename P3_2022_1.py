import numpy as np

def calcular_xn(x0, max_iter):
    """
    Calcular los siguientes valores Xn+1 utilizando la fórmula Xn+1 = 3.9 * Xn * (1 - Xn).

    Parametros:
        x0 (float): Valor inicial.
        max_iter (int): Número máximo de iteraciones.

    Returns:
        numpy.ndarray: Arreglo con los valores Xn+1 calculados.
    """
    xn_values = np.zeros(max_iter)
    xn_values[0] = x0

    for i in range(1, max_iter):
        xn_values[i] = 3.9 * xn_values[i-1] * (1 - xn_values[i-1])

    return xn_values

# Parámetros iniciales
x0 = 0.5  # Valor inicial
max_iter = 21  # Número máximo de iteraciones

# Calcular Xn+1
xn_values = calcular_xn(x0, max_iter)

print('n\tX_n')

# Mostrar los valores calculados
for i, xn_plus_1 in enumerate(xn_values):
    print(f'{i}: {xn_plus_1:.10f}')
