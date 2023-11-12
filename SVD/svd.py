import numpy as np


def calcMat(M, opc):
    if opc == 1:
        newM = np.dot(M.T, M)
    if opc == 2:
        newM = np.dot(M, M.T)

    eigenvalues, eigenvectors = np.linalg.eig(newM)
    ncols = np.argsort(eigenvalues)[::-1]

    if opc == 1:
        return eigenvectors[:, ncols].T

    else:
        return eigenvectors[:, ncols]


def calcD(M):
    if (np.size(np.dot(M, M.T))) > np.size(np.dot(M.T, M)):
        newM = np.dot(M.T, M)
    else:
        newM = np.dot(M, M.T)

    eigenvalues, eigenvectors = np.linalg.eig(newM)
    eigenvalues = np.sqrt(eigenvalues)

    return eigenvalues[::-1]


A = np.array([[0.92, 0.08], [0.08, 0.92]])
Vt = calcMat(A, 1)
U = calcMat(A, 2)
Sigma = calcD(A)
print("VT:")
print(Vt, "\n")
print("U")
print(U, "\n")
print("D")
print(Sigma)
