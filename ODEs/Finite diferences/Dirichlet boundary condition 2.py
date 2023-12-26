import numpy as np
import matplotlib.pyplot as plt
from Ej1 import soltrinfcol, soltrsupfil

def two_points(f, a, b, Ua, Ub, n, x0, x1):

    A = np.eye(n-2) * (-2)
    z = np.ones(n-3)
    upper = np.diag(z,1)
    lower = np.diag(z,-1)
    A += upper + lower

    h = 1/n
    A = A * (- 1/(h**2))  # El menos es para hacer definida positiva a la matriz
    
    F = np.linspace(a + 1/(n-2), b - 1/(n-2), n-2)

    k = 0
    for i in F:
        F[k] = f(i, x0, x1)
        k += 1

    F[0] = F[0] - Ua/(h**2)
    F[-1] = F[-1] - Ub/(h**2)

    G = np.linalg.cholesky(A)
    y = soltrinfcol(G,F)
    u = soltrsupfil(G.T,y)
    U = np.zeros(n)
    U[0], U[-1] = Ua, Ub
    U[1:n-1] = - u # Multiplico por menos para obtener la solucion correcta
    return U

#a
'''
fun = lambda x: - (np.pi**2) * np.cos(np.pi * x)
u = two_points(fun, 0, 1, 0, -1, 50)
print(u)
fig, ax = plt.subplots(1,1)
x = np.linspace(0,1,50)
ax.plot(x, u, color='purple', label='u(x)')
ax.grid()
ax.legend()
plt.show()
'''

#b

def f(x,x0,x1):
    if x > x0 and x < x1:
        return 1
    else:
        return 0


tuplas = [[0.25, 0.75], [0.1, 0.9], [0.2,0.5], [0.8,0.95], [0.05,0.95], [0.48, 0.52]]
x = np.linspace(0,1,50)
fig, ax = plt.subplots(1,1)
for i in tuplas:
    u = two_points(f, 0, 1, 0, 0, 50, i[0], i[1])
    ax.plot(x, u, label=f'Solucion con {i}')
ax.grid()
ax.legend()
plt.show()