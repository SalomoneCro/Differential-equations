import numpy as np
'''
def Runge_Kutta_Sistema(v_fun, interval, X0, h):
    X_i = X0
    t_i = interval[0]
    its =  int(np.abs(interval[0] - interval[1]) / np.abs(h))
    m = len(X0)
    sol = np.zeros((m,its))
    sol[:,0] = X0

    for i in range(its - 1):
        F1 = h*v_fun(t_i, X_i)
        F2 = h*v_fun(t_i + h/2, X_i + F1/2)
        F3 = h*v_fun(t_i + h/2, X_i + F2/2)
        F4 = h*v_fun(t_i + h, X_i + F3)
        X_i += (F1 + 2*F2 + 2*F3 + F4)/6
        t_i += h
        sol[:, i + 1] = X_i

    return sol
'''
def Runge_Kutta_Sistema(fun, t0, tf, X0, h):  #ahora es una lista!
    
    X0=np.array(X0)
    aproximaciones = np.array(X0)
    lista_tiempo = [t0]
    
    total_pasos = int((tf-t0)/h) 

    for i in range(total_pasos):
        F1 = h * fun(t0,X0)
        F2 = h * fun(t0 +.5*h, X0 + .5*F1)
        F3 = h * fun(t0 +.5*h, X0 + .5*F2)
        F4 = h * fun(t0 + h, X0 + F3)
        
        X_nuevo = X0 + (1/6)*( F1 + 2*F2 + 2*F3 + F4 ) 
        t0 = t0+h
        X0=X_nuevo  
        
        aproximaciones = np.vstack([aproximaciones,X0])
        lista_tiempo.append(t0)

    return aproximaciones

def Lines(a, b, T, U0, g1, g2, hx, ht):

    n = int((b-a) / hx)

    A = np.eye(n) * (-2)
    z = np.ones(n-1)
    upper = np.diag(z,1)
    lower = np.diag(z,-1)
    A += upper + lower
    A = A / (hx**2)
    
    def Ut(t,X):
        return -np.sin(t) * X**2 * np.sin(np.pi*X)
    
    def Uxx(t,X):
        return np.cos(t)*(2*np.sin(np.pi*X) + 4*np.pi*X*np.cos(np.pi*X) - (np.pi*X)**2 *np.sin(np.pi*X))
    
    def f(t, X):
        return -Uxx(t, X) + Ut(t,X)

    def g(t):
        v = np.zeros(n)
        v[0] = g1(t)
        v[-1] = g2(t)
        return v / (hx**2)
    
    grilla_espacial = np.arange(a,b,hx)
    k = 0
    cond_inicial = np.zeros(n)
    for i in grilla_espacial:
        cond_inicial[k] = U0(i)
        k += 1

    def F(t, X):
        return np.array(A @ X + g(t) + f(t, grilla_espacial))
    
    sol = Runge_Kutta_Sistema(F, 0, T, cond_inicial, ht)

    return sol


def Crank_Nicolson(a, b, T, U0, g1, g2, hx, ht):

    n = int((b-a) / hx)
    m = int(T/ht)

    alfa = ht/(2*(hx**2))

    A = np.eye(n) * (2*alfa+1)
    z = np.ones(n-1)
    upper = np.diag(z,1) * (-alfa)
    lower = np.diag(z,-1) * (-alfa)
    A += upper + lower
    B = np.eye(n) * (-2*alfa+1)
    B -= (upper + lower)


    def g(t):
        v = np.zeros(n)
        v[0] = (g1(t) + g1(t+ht)) * alfa
        v[-1] = (g2(t) + g2(t+ht)) * alfa
        return v 
    
    def Ut(t,X):
        return -np.sin(t) * X**2 * np.sin(np.pi*X)
    
    def Uxx(t,X):
        return np.cos(t)*(2*np.sin(np.pi*X) + 4*np.pi*X*np.cos(np.pi*X) - (np.pi*X)**2 *np.sin(np.pi*X))
    
    def f(t, X):
        return -Uxx(t, X) + Ut(t,X)
    
    def F(t,X):
        return ht * (f(t,X) + f(t+ht,X)) / 2 + g(t)
    
    
    grilla_espacial = np.arange(a,b,hx)
    sol = np.zeros((m + 1, n))

    k = 0
    for i in grilla_espacial:
        sol[0,k] = U0(i)
        k += 1

    t = ht
    for i in range(1, m+1):
        
        sol[i,:] = np.linalg.solve(A, B @ sol[i-1, :] + F(t,grilla_espacial))
        t += ht

    return sol

U0 = lambda x : x**2 * np.sin(np.pi*x)
g1 = lambda t : 0
g2 = lambda t : 0

a = 0
b = 1
T = 1
hx = 0.025
ht = 0.025
'''
a = 0
b = 1
T = 1
hx = 0.05
ht = 0.025*hx
'''
#sol = Lines(a,b,T,U0,g1,g2,hx,ht)
sol = Crank_Nicolson(a,b,T,U0,g1,g2,hx,ht)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
sol_an = lambda x,t: np.cos(t) * x**2 * np.sin(x*np.pi)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # Use add_subplot to create 3D axes

y = np.linspace(0, T, int(T/ht)+1)
x = np.linspace(a, b, int((b-a)/hx))
X, Y = np.meshgrid(x, y)  # grid of points

surf = ax.plot_surface(X, Y, sol, rstride=1, cstride=1,
                      cmap=cm.RdBu, linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y, sol_an(X, Y), rstride=1, cstride=1,
                      cmap=cm.RdBu, linewidth=0, antialiased=False)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()



'''
import pylab
x = np.arange(-3.0,3.0,0.1)
y = np.arange(-3.0,3.0,0.1)
X,Y = pylab.meshgrid(x, y) # grid of point
Z = sol_an(X, Y) # evaluation of the function on the grid
im = pylab.imshow(Z,cmap=pylab.cm.RdBu) # drawing the function
# adding the Contour lines with labels
cset = pylab.contour(Z,np.arange(-1,1.5,0.2),linewidths=2,cmap=pylab.cm.Set2)
pylab.clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
pylab.colorbar(im) # adding the colobar on the right

pylab.show()
'''