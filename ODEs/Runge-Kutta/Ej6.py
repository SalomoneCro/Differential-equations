import numpy as np
import matplotlib.pyplot as plt 
import sympy as sp
from math import factorial

def Runge_Kutta(fun, interval, x_0, h):
    x_i = x_0
    t_i = interval[0]
    its =  int(np.abs(interval[0] - interval[1]) / np.abs(h))
    sol = [x_0]

    for _ in range(its):

        F1 = h*fun(t_i, x_i)
        F2 = h*fun(t_i + h/2, x_i + F1/2)
        F3 = h*fun(t_i + h/2, x_i + F2/2)
        F4 = h*fun(t_i + h, x_i + F3)
        x_i += (F1 + 2*F2 + 2*F3 + F4)/6
        t_i += h
        sol.append(x_i)
    return sol

fun = lambda t, x: - 5 * x 
t = sp.Symbol('t')
x = sp.Function('x')(t)
funs = - 5 * x 

def taylor_method(order, f, x_0, h, interval):
    t = sp.symbols('t')
    x = sp.Function ("x")(t)
    x_sol = [x_0]
    func = f
    f_der=[f] #List of derivatives

    for i in range(order-1):  #Calculates f', f'', ...
        f = sp.diff(f,t)
        f = f.xreplace({sp.diff(x,t):func}) #replace dx/dt with f
        f_der.append(f)

    #steps =[a, a+h, a+2h, ..., b]
    a = interval[0]
    b = interval[1]
    steps = np.arange(a,b,h,dtype=float) if a!=b else [a]
    np.append(steps, b)

    for j in steps[1:]: #Iterates for all points
        x_1 = x_0
        for n in range(order): #calculates x(t + h)
            der_value = f_der[n].xreplace({x: x_1, t:j})
            x_1 += (h**(n+1)) * der_value / factorial(n+1)
        x_0 = x_1
        x_sol.append(x_0)
    return x_sol, steps

step = 0.2

solsEuler = taylor_method(1, funs, 1, step, [0,1])[0]
solsRK = Runge_Kutta(fun, [0,1], 1, step)  

space = np.linspace(0,1,int(1/step) + 1)

fig, ax = plt.subplots(1,1)

ax.plot(space, solsRK, label='Solucion R-K')
#ax.plot(space, solsEuler, label='Solucion Euler')
ax.plot(space, np.exp(-5 * space), label='Solucion Analitica')
ax.legend()
ax.grid()
plt.show()