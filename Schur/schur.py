import numpy as np


def schur_decomposition(A, max_iterations=100):
    n = A.shape[0]
    Q = np.eye(n)
    T = A.copy()

    for _ in range(max_iterations):
        Q_, R = np.linalg.qr(T)

        T = np.dot(R, Q_)

        Q = np.dot(Q, Q_)

        if np.allclose(T, np.triu(T)):
            break

    return Q, T


A = np.array([[2, -1], [1, 3]])

Q, T = schur_decomposition(A)

print("Matriz triangular superior (T):")
print(T)

print("Matriz unitaria (Q):")
print(Q)
