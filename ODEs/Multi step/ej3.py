import numpy as np
import matplotlib.pyplot as plt

def Multipaso(fun, interval, x_0, h):
    x_i = x_0
    t_is = [interval[0]]
    its =  int(np.abs(interval[0] - interval[1]) / np.abs(h))
    sol = [x_0]

    for _ in range(3):

        F1 = h*fun(t_is[-1], x_i)
        F2 = h*fun(t_is[-1] + h/2, x_i + F1/2)
        F3 = h*fun(t_is[-1] + h/2, x_i + F2/2)
        F4 = h*fun(t_is[-1] + h, x_i + F3)
        x_i += (F1 + 2*F2 + 2*F3 + F4)/6
        sol.append(x_i)
        t_is.append(t_is[-1] + h)

    sol.append(sol[-1] + h/24*(55*fun(t_is[-1], sol[-1]) - 59*fun(t_is[-2], sol[-2]) + 
                                37*fun(t_is[-3], sol[-3]) - 9*fun(t_is[-4], sol[-4])))
                                
    return sol


fun = lambda t, x: -2*t*(x**2)

solRKM = Multipaso(fun, [0,1], 1, 0.25)
solan = lambda x: 1/(1+x**2)

x = np.linspace(0,1,5)
fig, ax = plt.subplots(1,1)
ax.plot(x, solan(x), label='Analitica')
ax.plot(x, solRKM, label='Runge-Kutta Multipaso')
ax.grid()
ax.legend()
plt.show()