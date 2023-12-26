import numpy as np
import matplotlib.pyplot as plt

fun1 = lambda x1, x2, t: x2
fun2 = lambda x1, x2, t: -192 * x1

def F(t, X):
    return np.array([fun1(X[0], X[1], t), fun2(X[0], X[1], t)])

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

step = 0.01 
sols = Runge_Kutta_Sistema(F, [0,5], [1/6,0], step)



fig, ax = plt.subplots(1,1)
space = np.linspace(0,5,int(5/step))
ax.plot(space, sols[:,0], label='Solucion', color='purple')
ax.grid()
ax.legend()
plt.show()

