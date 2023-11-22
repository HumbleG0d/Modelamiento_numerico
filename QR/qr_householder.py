import numpy as np


def column_convertor(x):
    """
    Converts 1d array to column vector
    """
    x.shape = (1, x.shape[0])
    return x


def get_norm(x):
    """
    Returns Norm of vector x
    """
    return np.sqrt(np.sum(np.square(x)))


def householder_transformation(v, iter):
    """
    Returns Householder matrix for vector v
    """
    size_of_v = v.shape[1]
    e1 = np.zeros_like(v)
    e1[0, 0] = 1
    vector = get_norm(v) * e1
    if v[0, 0] < 0:
        vector = -vector
    u = (v + vector).astype(np.float32)
    H = np.identity(size_of_v) - (
        (2 * np.matmul(np.transpose(u), u)) / np.matmul(u, np.transpose(u))
    )
    print(f"H{iter+1}: ", H)
    return H


def qr_step_factorization(q, r, iter, n):
    """
    Return Q and R matrices for iter number of iterations.
    """
    v = column_convertor(r[iter:, iter])
    Hbar = householder_transformation(v, iter)
    H = np.identity(n)

    H[iter:, iter:] = Hbar
    r = np.matmul(H, r)
    q = np.matmul(q, H)
    return q, r


def QR_TH(A):
    n, m = A.shape
    Q = np.identity(n)
    R = A.astype(np.float32)
    for i in range(min(n, m)):
        # For each iteration, H matrix is calculated for (i+1)th row
        Q, R = qr_step_factorization(Q, R, i, n)
    min_dim = min(m, n)
    R = np.around(R, decimals=6)
    R = R[:min_dim, :min_dim]
    Q = np.around(Q, decimals=6)

    return Q, R


A = np.array([[1, 1, 9], [1, 1, 1], [2, -1, 0]])

b = np.array([[436], [44], [0]])

print("Planteamiento del sistema")
print("Matriz A: ")
print(A)

print("Matriz b: ")
print(b)

print("===================================")

Q, R = QR_TH(A)

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
