import numpy as np
import matplotlib.pyplot as plt

def ghost_at_b(f, a, b, Ua, Uxb, n):

    h = (b - a) / n

    #Defino la grilla
    gri = np.linspace(a,b,n) 

    #Defino la matriz
    A = np.zeros((n - 1, n - 1))  # No es n-2 porque tengo una condicion de Neumann, debo despejar el valor en un extremo
    
    A[0,0] = -2/(h**2)
    A[0,1] = 1/(h**2)
    
    for i in range(1, n-2):
        A[i, i-1] = A[i-1, i]
        A[i,i] = -2/(h**2)
        A[i,i+1] = 1/(h**2)

    #Agrego la ultima fila
    A[n-2, n-3] = 1/(h**2)
    A[n-2, n-2] = -3/(h**2)

    #Defino el vector b de Au = b

    b = np.zeros(n-1)

    b[0] = f(gri[1]) - Ua / (h**2)
    for i in range(1, n-2):
        b[i] = f(gri[i])
    b[n-2] = f(gri[n-2]) - 2*Uxb/h

    U = np.linalg.solve(A,b)
    u = np.zeros(n)
    u[0] = Ua
    u[1:] = U
    return u, gri
