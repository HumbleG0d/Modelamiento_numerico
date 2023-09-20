import numpy as np

def calcular(x0 , max_iter):
  """
    Tenemos que x0 = ln(1.2)
    x_k = 1/k - 5x_k-1
    
    Parametros:
        x0 (float): Valor inicial.
        max_iter (int): Número máximo de iteraciones.

  """

  xk_values = np.zeros(max_iter)
  xk_values[0] = x0

  for i in range(1 , max_iter):
    xk_values[i] = 1/i - 5*xk_values[i-1]
  
  return xk_values

#Valor inicial
x0 = np.log(1.2)

#Nùmero de iteraciones
max_iter = 20

#Calculamos los valores
xk_values = calcular(x0 , max_iter)

print('k\tx_k')
for i , xk in enumerate(xk_values , start = 1):
  print(f'{i}: {xk: .10f}')