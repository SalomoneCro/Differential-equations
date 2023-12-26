import sympy as sp
import numpy as np
from math import factorial

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