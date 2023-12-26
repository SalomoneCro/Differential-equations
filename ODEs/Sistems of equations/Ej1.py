import numpy as np
import matplotlib.pyplot as plt


def F(t, x):
    return np.array([t + x[0]**2 + x[1], t**2 - x[0] + x[1]**2],dtype=float)

def Runge_Kutta_Sistema(v_fun, interval, X0, h):
    X_i = X0
    t_i = interval[0]
    its =  int(np.abs(interval[0] - interval[1]) / np.abs(h))
    sol = [X0]

    for _ in range(its - 1):
        F1 = h*v_fun(t_i, X_i)
        F2 = h*v_fun(t_i + h/2, X_i + F1/2)
        F3 = h*v_fun(t_i + h/2, X_i + F2/2)
        F4 = h*v_fun(t_i + h, X_i + F3)
        X_i += (F1 + 2*F2 + 2*F3 + F4)/6
        t_i += h
        sol.append(X_i.tolist())
    return np.array(sol)

sols = Runge_Kutta_Sistema(F, [-1,1], [0.43,-0.69], 0.01)

fig, ax = plt.subplots(1,1)
x = np.linspace(1,2,200)
ax.plot(x, sols[:,0])
ax.plot(x, sols[:,1])
ax.grid()
plt.show()
