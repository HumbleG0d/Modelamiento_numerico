import numpy as np


def gram_schmidt_qr(A):
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for i in range(n):
        v = A[:, i]
        for j in range(i):
            R[j, i] = np.dot(Q[:, j], A[:, i])
            print(f"R[{j},{i}] : ", R[j, i])
            v = v - R[j, i] * Q[:, j]
            print("v:", v)
        R[i, i] = np.linalg.norm(v)
        print(f"R[{i},{i}] : ", R[i, i])
        Q[:, i] = v / R[i, i]
        print(f"Q{i} : ", Q[:, i])

    return Q, R


# Ejemplo de uso:
A = np.array(
    [
        [4, -1, 0, -1, 0, 0],
        [-1, 4, -1, 0, -1, 0],
        [0, 1, 4, 0, 0, -1],
        [-1, 0, 0, 4, -1, 0],
        [0, -1, 0, -1, 4, -1],
        [0, 0, -1, 0, -1, 4],
    ]
)
b = np.array([[56], [25], [62], [41], [10], [47]])
print("Planteamiento del sistema")
print("Matriz A: ")
print(A)

print("Matriz b: ")
print(b)

print("===================================")

Q, R = gram_schmidt_qr(A)

print("\nMatriz Q (ortogonal):")
print(Q)
print("\nMatriz R (triangular superior):")
print(R)

print(
    """
      Tenemos el sistema Ax=b.
      Mediante cholesKy tenemos que A=RQ
      Entonces nuestro nuevo sistema es RQx = b
      Para resolver dicho sistema Rx=Q^Tb
      sea Q^Tb = Y
"""
)

Y = np.dot(Q.T, b)
print("Matriz Y : ")
print(Y)

n = A.shape[1]
print("\nCalculando la matriz solucion(x) por sustitucion regresiva.")
x = np.zeros(n)

for i in range(n - 1, -1, -1):
    x[i] = (Y[i] - np.dot(R[i, i + 1 :], x[i + 1 :])) / R[i, i]
print("Matriz solucion:")
print(x)
