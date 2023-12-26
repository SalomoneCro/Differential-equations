import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return -1

def dif_centrada(f, a, b, Ua, Ub, n, ep):

    h = (b - a) / n

    #Defino la grilla
    gri = np.linspace(a,b,n) 

    #Defino la matriz
    A = np.zeros((n - 2, n - 2))
    
    A[0,0] = -2*ep/(h**2)
    A[0,1] = (2*ep - h) / (2*(h**2))

    A[-1,-2] = (2*ep + h) / (2*(h**2))
    A[-1,-1] = -2*ep/(h**2)
    
    for i in range(1, n-3):
        A[i, i-1] = (2*ep + h) / (2*(h**2))
        A[i,i] = -2*ep/(h**2)
        A[i,i+1] = (2*ep - h) / (2*(h**2))

    #Defino el vector b de Au = b

    b = np.zeros(n-2)

    b[0] = f(gri[1]) - (2*ep+h) / (2*(h**2))
    for i in range(1, n-2):
        b[i] = f(gri[i])
    b[n-3] = f(gri[-1]) - (6*ep-3*h) / (2*(h**2))

    U = np.linalg.solve(A,b)
    u = np.zeros(n)
    u[0] = Ua
    u[1:-1] = U
    u[-1] = Ub
    return u, gri

def upwind(f, a, b, Ua, Ub, n, ep):

    h = (b - a) / n

    #Defino la grilla
    gri = np.linspace(a,b,n) 

    #Defino la matriz
    A = np.zeros((n - 2, n - 2))
    
    A[0,0] = -(2*ep + h)/(h**2)
    A[0,1] = ep/(h**2)

    A[-1,-2] = (ep + h) / (h**2)
    A[-1,-1] = -(2*ep + h)/(h**2)
    
    for i in range(1, n-3):
        A[i, i-1] = (ep + h) / (h**2)
        A[i,i] = -(2*ep + h)/(h**2)
        A[i,i+1] = ep / (h**2)

    #Defino el vector b de Au = b

    b = np.zeros(n-2)

    b[0] = f(gri[1]) - (ep+h) / (h**2)
    for i in range(2, n-2):
        b[i] = f(gri[i])
    b[n-3] = f(gri[-1]) - 3*ep / (h**2)

    U = np.linalg.solve(A,b)
    u = np.zeros(n)
    u[0] = Ua
    u[1:-1] = U
    u[-1] = Ub
    return u, gri

eps = [0.3, 0.1, 0.05, 0.0005]

fig, ax = plt.subplots(2,1)

for i in eps:
    U, x = dif_centrada(f, 0, 1, 1, 3, 100, i)
    Uu, xx = upwind(f, 0, 1, 1, 3, 80, i)
    ax[0].plot(x, U, label=f'DFC eps={i}')
    ax[1].plot(xx, Uu, label=f'upwind eps={i}')
ax[0].grid()
ax[0].legend()
ax[1].grid()
ax[1].legend()
plt.show()