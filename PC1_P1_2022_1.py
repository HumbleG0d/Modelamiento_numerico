import numpy as np

def calcular(max_iter):
  """
    x_n = sen(1/n)/(1/n) ; n>=1  
  """

  x_values = np.zeros(max_iter)

  for i in range(1, max_iter):
    x_values[i] = np.sin(1/i)/(1/i)
  
  return x_values
#Nùmero de iteraciones
max_iter = 11

x_values = calcular(max_iter)
print(x_values)

print("n\tx_n")

for i , xn in enumerate(x_values):
  if(i != 0):  print(f'{i} : {xn : .10f}')

