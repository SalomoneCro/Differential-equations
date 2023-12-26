import numpy as np
import matplotlib.pyplot as plt
import time

def soltrinfcol(A,h):
    b = h.copy()
    x = np.zeros(len(h))
    k = np.min(np.nonzero(b))
    for j in range(k, len(b)):
        x[j] = b[j] / A[j,j]
        b[j:] = b[j:] - x[j] * A[j: ,j] 
    return x


def soltrsupfil(Z,c):
    b = c.copy()
    k = np.max(np.nonzero(b))
    for i in reversed(range(k + 1)):
        b[i] = (b[i] - Z[i, i + 1:] @ b[i + 1:] ) / Z[i,i]
    return b

def two_points(f, a, b, Ua, Ub, n):

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
        F[k] = f(i)
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


#Prueba de velocidad de algoritmo solo con solve de numpy
'''
def two_points2(f, a, b, Ua, Ub, n):

    A = np.eye(n-2) * (-2)
    z = np.ones(n-3)
    upper = np.diag(z,1)
    lower = np.diag(z,-1)
    A += upper + lower

    h = 1/n
    A = A * (1/(h**2))  
    
    F = np.linspace(a + 1/(n-2), b - 1/(n-2), n-2)

    k = 0
    for i in F:
        F[k] = f(i)
        k += 1

    F[0] = F[0] - Ua/(h**2)
    F[-1] = F[-1] - Ub/(h**2)

    u = np.linalg.solve(A,F)
    U = np.zeros(n)
    U[0], U[-1] = Ua, Ub
    U[1:n-1] = u
    return U


fun = lambda x: - (np.pi**2) * np.cos(np.pi * x)
x = time.time()
u = two_points(fun, 0, 1, 0, -1, 5000)
print(time.time() - x)
x = time.time()
u = two_points2(fun, 0, 1, 0, -1, 5000)
print(time.time() - x)
'''

