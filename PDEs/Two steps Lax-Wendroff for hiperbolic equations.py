import numpy as np

def LW2p(F, U0, Interval, T, ht, hx, g1):
    n = int(np.abs(Interval[1] - Interval[0])/hx)
    m = int(T/ht)

    grilla = np.linspace(Interval[0], Interval[1], n)
    def Ujunmedio(ui_act, ui_post):
        return (ui_post+ui_act)/2 - ht/(2*hx) * (F(ui_post) - F(ui_act))
    
    def Uj_unmedio(ui_ant, ui_act):
        return (ui_act+ui_ant)/2 - ht/(2*hx) * (F(ui_act) - F(ui_ant))

    sol = np.zeros((m,n))
    
    sol[0,:] = U0(grilla,n)
    for j in range(m): sol[j,0] = g1(j)
    for j in range(m-1):
        for i in range(1,n-1):
            sol[j+1,i] = sol[j,i] - ht/hx * (F(Ujunmedio(sol[j,i], sol[j,i+1])) - F(Uj_unmedio(sol[j,i-1],sol[j,i])))
        sol[j+1,n-1] = sol[j,n-1] = sol[j,n-1] - ht/hx * (F(Ujunmedio(sol[j,n-1], sol[j,n-1])) - F(Uj_unmedio(sol[j,n-2],sol[j,n-1])))
   
    return sol

def F(x): return (x**2)/2
a = -1
b = 1
I = [a,b]
T = 1
ht = 0.01 
hx = 0.08
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

g1 = lambda T:0


sol = LW2p(F,U01,I,T,ht,hx,g1)

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
