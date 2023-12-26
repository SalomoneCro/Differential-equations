import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

def diferencias_finitas(f,x):
    n = len(x)
    h = np.abs(x[-1] - x[0])/n

    #Defino la matriz de rigidez
    A = np.zeros((n-2,n-2))
    A = A + np.eye(n-2) * 2 / h + np.diag(np.ones(n-3),1) * (-1/h) + np.diag(np.ones(n-3),-1) * (-1/h)
    A = A + np.diag(np.ones(n-3),1) * (h/6) + np.diag(np.ones(n-3),-1) * (h/6)

    diag = np.eye(n-2)
    for i in range(1,n-1):
        diag[i-1,i-1] *= 2*h/3

    A = A + diag

    #Defino las funciones sombrero
    def phi(x,i,h,grid):
        if grid[i-1] <= x <= grid[i]:
            return (x - grid[i-1]) / h
        elif grid[i] <= x <= grid[i+1]:
            return (grid[i+1]-x) / h
        else:
            return 0

    #Esta función me definira cada entrada del vector de carga
    def fun(x,i,h,grid): 
        return f(x)*phi(x,i,h,grid)

    #Creo el vector de carga
    F=np.array([sc.integrate.quad(fun,x[i-1],x[i+1] ,args=(i,h,x))[0] for i in range(1,n-1)])

    sol = np.zeros(n)
    sol[1:n-1] = np.linalg.solve(A,F)
    return sol

#Parametros
fun = lambda x: np.sin(np.pi*x) * (1+np.pi**2)
sol_an = lambda x: np.sin(np.pi*x)
x = np.linspace(0,1,10)
sol = diferencias_finitas(fun, x)

#Grafico
ax, fig = plt.subplots(1,1)
plt.plot(x, sol, color='darkgreen', label='Sol con EF')
plt.plot(x, sol_an(x), color='purple', label='Sol exacta')
plt.legend()
plt.grid()
plt.show()