import numpy as np
import matplotlib.pyplot as plt 

def Runge_Kutta(fun, interval, x_0, h):
    x_i = x_0
    t_i = interval[0]
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

fun = lambda t, x: np.exp(x*t) + np.cos(x - t)

interMiami = [1, 1.04]

sols = Runge_Kutta(fun, interMiami, 3, 0.001)

space = np.linspace(1,1.04,40)

fig, ax = plt.subplots(1,1)

ax.plot(space, sols, label='Solucion R-K')
ax.legend()
ax.grid()
plt.show()