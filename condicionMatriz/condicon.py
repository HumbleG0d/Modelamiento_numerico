# Calcular la condicion de una matriz
import numpy as np

matriz = np.array(
    [
        [625, 125, 25, 5, 1],
        [256, 64, 16, 4, 1],
        [81, 27, 9, 3, 1],
        [16, 8, 4, 2, 1],
        [1, 1, 1, 1, 1],
    ]
)

try:
    matriz_inversa = np.linalg.inv(matriz)
    print("Matriz inversa:")
    print(matriz_inversa)
except np.linalg.LinAlgError:
    print("La matriz no es invertible.")

suma_maxima_filas = np.max(np.sum(np.abs(matriz), axis=1))
print("Norma infinita de la matriz A: ", suma_maxima_filas)
suma_maxima_filas_inv = np.max(np.sum(np.abs(matriz_inversa), axis=1))
print("Norma infinita de la inversa de la matriz A: ", suma_maxima_filas_inv)

numero_condicion = suma_maxima_filas * suma_maxima_filas_inv
print("Condicion :", numero_condicion)

# Calculo del Er
"""
  NornmaInfinita(R)/(Cond(A)*NormaInfinita(b)) <= Er <= NornmaInfinita(R)*Cond(A)/NormaInfinita(b)
  R = Ax - b
"""
