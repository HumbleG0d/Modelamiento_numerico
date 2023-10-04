def eliminacionGaussiana(A, b):
    x = [0] * len(b)
    M = A.copy()

    for i in range(len(b)):
        M[i].insert(len(b) + 1, b[i][0])

    print("Matriz Aumentada:\n")
    print("\n".join(["".join(["{:4}".format(item) for item in row]) for row in M]))
    print("\n")

    i = 0
    aux = False

    while i < len(A):
        if M[i][i] == 0:
            aux = True
            break
        else:
            for k in range(i + 1, len(A)):
                num = -M[k][i]
                for j in range(len(b) + 1):
                    M[k][j] = M[k][j] + num * (M[i][j] / M[i][i])

                print(
                    "\n".join(
                        [
                            "".join(["{:10.3f}".format(item) for item in row])
                            for row in M
                        ]
                    )
                )
                print("\n")
        i = i + 1

    print("\nMatriz triangular")
    print("\n".join(["".join(["{:10.3f}".format(item) for item in row]) for row in M]))

    for i in range(len(A)):
        for j in range(len(A[0])):
            A[i][j] = M[i][j]
    for i in range(len(b)):
        b[i][0] = M[i][len(b)]

    if aux:
        print("No hay solución")
    else:
        sustitucionRegresiva(A, b)


def sustitucionRegresiva(U, b):
    x = [0] * len(b)
    i = len(b) - 1
    num = 0

    while i >= 0:
        j = i + 1
        sum = 0
        if j != len(U):
            while j < len(U):
                sum = sum + x[j] * U[i][j]
                j = j + 1
        num = (b[i][0] - sum) / U[i][i]
        x[i] = num
        i = i - 1

    print("\nSolucion: ")
    x_round = [round(i, 2) for i in x]
    print(x_round)


A = [
    [1, 1, 1, 1],
    [120, 200, 80, 120],
    [-1, 1, 0, -1],
    [0, 1, -3, 0],
]
b = [
    [14],
    [2020],
    [0],
    [0],
]

eliminacionGaussiana(A, b)
