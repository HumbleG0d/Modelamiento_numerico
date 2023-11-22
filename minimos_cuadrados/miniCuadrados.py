import numpy as np


def are_columns_linearly_independent(matrix):
    matrix = matrix.astype(float)  # Convertir a tipo de datos flotante
    original_rank = np.linalg.matrix_rank(matrix)
    cleaned_rank = np.linalg.matrix_rank(
        np.unique(matrix, axis=1), tol=np.finfo(matrix.dtype).eps * 10
    )

    return original_rank == cleaned_rank


# Ejemplo de uso:
A = np.array([[1, -6], [1, -2], [1, 1], [1, 7]])
b = np.array([[-1], [2], [1], [6]])
value = are_columns_linearly_independent(A)
print("¿Las columnas de matrix1 son linealmente independientes?", value)
if value:
    print("Entonces A^TA es inversible")
    """
    Tenemos el siguiente sistema A^TAx = A^Tb
    como A^TA es inversible , entonces x = [(A^TA)^-1]A^Tb
    Sean : Q = (A^TA)^-1 y W = A^Tb
    Asi x = QW
  """
    Q = np.linalg.inv(np.dot(A.T, A))
    print(Q)
    W = np.dot(A.T, b)
    print(W)
    x = np.dot(Q, W)
    print("Matriz solución: ")
    print(x)
else:
    print("El sistema no se puede resolver por minimos cuadrados")
