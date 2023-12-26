import numpy as np
import matplotlib.pyplot as plt

fun1 = lambda x1, x2, t: x1**(-2) + np.log(x2) + t**2
fun2 = lambda x1, x2, t: np.exp(x2) - np.cos(x1) + np.sin(t)*x1 - (x1*x2)**(-3)

def F(t, X):
    return np.array([fun1(X[0], X[1], t), fun2(X[0], X[1], t)])

def Runge_Kutta_Sistema(v_fun, interval, X0, h):
    X_i = X0
    t_i = interval[-1]
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

sols = Runge_Kutta_Sistema(F, [1,2], [-2,1], -0.01)
print(sols)
fig, ax = plt.subplots(1,1)
x = np.linspace(1,2,100)
ax.plot(x, np.flip(sols[:,0]))
ax.plot(x, np.flip(sols[:,1]))
ax.grid()
plt.show()

