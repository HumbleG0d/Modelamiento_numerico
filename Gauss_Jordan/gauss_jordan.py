import numpy as np
import sys

def gauss_jordan(A,b):
    M = A.copy()
    n = len(b)
    x = np.zeros(n)

    for i in range(len(b)):
        M[i].insert(len(b) + 1 , b[i][0])
    
    print("Matriz Aumentada: ")
    print('\n'.join([''.join(['{:6}'.format(item) for item in row]) for row in M]))

    for i in range(n):
        if M[i][i] == 0.0:
            sys.exit('Division por 0 detectada')

        for j in range(n):
            if i != j:
                ratio = M[j][i] / M[i][i]

                for k in range(n+1):
                    M[j][k] = M[j][k] - ratio * M[i][k]

    print("\nMatriz Triangular:")
    print('\n'.join([''.join(['{:10.4f}'.format(item) for item in row]) for row in M]))

    for i in range(n):
        x[i] = M[i][n] / M[i][i]
    print("\nMatriz solución: ")
    for i in range(n):
        print('X%d = %0.2f' % (i , x[i]) , end = '\t')

A = [[1.5 , 2 ,3] ,
     [1 , 1.5 , 2],
     [1 , 1 ,1],]
b = [[250] ,
     [175],
     [140],]

gauss_jordan(A,b)