import numpy as np
import matplotlib.pyplot as plt 

def Runge_Kutta(fun, interval, x_0, h):
    x_i = x_0
    t_i = interval[1]
    its =  int(np.abs(interval[0] - interval[1]) / np.abs(h))
    sol = [x_0]

    for _ in range(its - 1):

        F1 = h*fun(t_i, x_i)
        F2 = h*fun(t_i + h/2, x_i + F1/2)
        F3 = h*fun(t_i + h/2, x_i + F2/2)
        F4 = h*fun(t_i + h, x_i + F3)
        x_i += (F1 + 2*F2 + 2*F3 + F4)/6
        t_i += h
        sol.append(x_i)
    return sol

fun = lambda t,x: (x - x*np.exp(t)) / (np.exp(t) + 1)

sols = Runge_Kutta(fun, [-2,0], 3, -0.01)


def funcion(t):
    return np.exp(t-2*np.log(np.exp(t)+1))*12

x = np.linspace(-2,0, 200)
fig, ax = plt.subplots(1,1)
ax.plot(x, np.flip(sols), label='Runge-Kutta')
ax.plot(x,funcion(x), label='Solucion analitica')
ax.grid()
ax.legend()
plt.show()


    