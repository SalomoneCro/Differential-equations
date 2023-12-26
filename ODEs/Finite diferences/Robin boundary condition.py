import numpy as np
import matplotlib.pyplot as plt
from math import e
def self_adjoint_bvp(f, p, q, a, b, Ua, s, t, d, n):
    h = (b - a) / n

    #Defino la grilla
    gri = np.linspace(a,b,n) 

    #Defino grilla 1 medios
    gri_medios = []
    for i in range(len(gri) - 1):
        gri_medios.append((gri[i] + gri[i+1])/2)

    #Defino la matriz
    A = np.zeros((n - 1, n - 1))  # No es n-2 porque tengo una condicion de robin, debo despejar el valor en un extremo
    
    A[0,0] = -(p(gri_medios[0]) + p(gri_medios[1])) / (h**2) - q(gri[1])
    A[0,1] = p(gri_medios[1]) / (h**2)
    
    for i in range(1, n-2):
        A[i, i-1] = A[i-1, i]
        A[i,i] = -(p(gri_medios[i]) + p(gri_medios[i+1])) / (h**2) - q(gri[i+1])
        A[i,i+1] = p(gri_medios[i+1]) / (h**2)

    #Agrego la ultima fila
    A[n-2, n-3] = -t/h
    A[n-2, n-2] = s + t/h

    #Defino el vector b de Au = b

    b = np.zeros(n-1)

    b[0] = f(gri[1]) - p(gri_medios[0]) * Ua / (h**2)
    for i in range(1, n-2):
        b[i] = f(gri[i])
    b[n-2] = d

    U = np.linalg.solve(A,b)
    u = np.zeros(n)
    u[0] = Ua
    u[1:] = U
    return u, gri

p = lambda x: 1 + x**2
q = lambda x: x
anali = lambda x: np.exp(-x) * ((x-1)**2)

def an(x):
    u=(e**(-x))*((x**2) -(2*x)+1)
    p1=2*x
    p=1 + x**2
    u1=(e**(-x))*(-(x**2) +(4*x)-3)
    u2=(e**(-x))*((x**2) -(6*x)+7)
    q=x
    return p1*u1 + p*u2 - q*u

U, x = self_adjoint_bvp(an,p,q,a=0,b=1,Ua=1,s=2,t=-3,d=0,n=80)

fig, ax = plt.subplots(1,1)
ax.plot(x, U, label='Solucion Numerica')
ax.plot(x, anali(x), label='Solucion Analitica')
ax.grid()
ax.legend()
plt.show()
