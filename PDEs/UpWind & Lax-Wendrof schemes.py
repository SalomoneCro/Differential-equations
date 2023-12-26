import numpy as np    
'''
def Upwind(U0, a, interval, T, ht, hx):
    n = int((interval[1] - interval[0]) / hx)
    m = int(T/ht)
    grilla_x = np.linspace(interval[0], interval[1], n)
    u = U0(grilla_x, n)

    beta = a*ht/hx
    gamma = 1 - beta

    A = np.eye(n) * (gamma)
    z = np.ones(n-1) * beta
    lower = np.diag(z,-1)
    A += lower

    sol = np.zeros((m, n))

    for i in range(m):
        sol[i,:] = A @ u
        u = sol[i,:]

    return sol

def Lax_Wendroff(U0, a, interval, T, ht, hx):
    n = int((interval[1] - interval[0]) / hx)
    m = int(T/ht)

    grilla_x = np.linspace(interval[0], interval[1], n)
    u = U0(grilla_x, n)

    c = a * ht / hx

    A = np.eye(n) * (1 - c**2)
    sub_diag = np.ones(n-1) * c * (1+c) / 2
    supra_diag = np.ones(n-1) * (c * (1-c) / 2)
    lower = np.diag(sub_diag,-1)
    upper = np.diag(supra_diag,1)
    A += lower - upper
    
    sol = np.zeros((m, n))
    sol[0,:] = u
    for i in range(1, m):
        sol[i,:-1] = (A @ u)[:-1]
        sol[i,-1] = c * (1+c) / 2 * u[-2] + (1 - c**2) * u[-1] - (c * (1-c) / 2) * u[-1]
        u = sol[i,:]

    return sol

a = -1
b = 1
T = 1
#ht tiene que ser menor que hx
ht = 0.02
hx = 0.03

def U01(x, l):
    lista = np.zeros(l)
    for i in range(len(x)):
        lista[i] = (x[i]+1)*np.exp(-x[i]/2)
    return lista

def U02(x, l):
    lista = np.zeros(l)
    for i in range(len(x)):
        if np.abs(x[i]) > 1/2:
            lista[i] = 0
        else:
            lista[i] = 1
    return lista

sol = Lax_Wendroff(U02, 1, [a,b], T, ht, hx)


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # Use add_subplot to create 3D axes

y = np.linspace(0, T, int(T/ht))
x = np.linspace(a, b, int((b-a)/hx))
X, Y = np.meshgrid(x, y)  # grid of points


surf = ax.plot_surface(X, Y, sol, rstride=1, cstride=1,
                      cmap=cm.RdBu, linewidth=0, antialiased=False)


ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
'''
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
def Lax_Wendroff(U0, a, interval, T, ht, hx):
    n = int((interval[1] - interval[0]) / hx)
    m = int(T/ht)

    grilla_x = np.linspace(interval[0], interval[1], n)
    u = U0(grilla_x, n)

    c = a * ht / hx

    A = np.eye(n) * (1 - c**2)
    sub_diag = np.ones(n-1) * c * (1+c) / 2
    supra_diag = np.ones(n-1) * (c * (1-c) / 2)
    lower = np.diag(sub_diag,-1)
    upper = np.diag(supra_diag,1)
    A += lower - upper
    
    sol = np.zeros((m, n))
    sol[0,:] = u
    for i in range(1, m):
        sol[i,:-1] = (A @ u)[:-1]
        sol[i,-1] = c * (1+c) / 2 * u[-2] + (1 - c**2) * u[-1] - (c * (1-c) / 2) * u[-1]
        u = sol[i,:]

    return sol

a = 0
b = 10
a0 = 1
T = 17
#ht tiene que ser menor que hx
ht = 0.05
hx = 0.08
def U(x, l):
    lista = np.zeros(l)
    for i in range(len(x)):
        lista[i] = np.exp(-20*(x[i]-2)**2) + np.exp(-(x[i]-5)**2)
    return lista

sol = Lax_Wendroff(U, a0, [a,b], T, ht, hx)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # Use add_subplot to create 3D axes

y = np.linspace(0, T, int(T/ht))
x = np.linspace(a, b, int((b-a)/hx))
X, Y = np.meshgrid(x, y)  # grid of points


surf = ax.plot_surface(X, Y, sol, rstride=1, cstride=1,
                      cmap=cm.RdBu, linewidth=0, antialiased=False)


ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
