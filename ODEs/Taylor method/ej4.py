import sympy as sp
from TaylorGral import taylor_method
import numpy as np

'''
t = sp.Symbol('t')
x = sp.Function('x')(t)

my_function = x + sp.exp(t) + t*x

def x1(x, ev0):
    return ev0 + exp(x) + x * ev0

    return ev1 + exp(x) + ev0 + x * ev1

def x3(x, ev1, ev2):
    return ev2 + exp(x) + 2*ev1 + x*ev2

def x4(x,ev2,ev3):
    return ev3 + exp(x) + 3*ev2 + x*ev3

def x5(x,ev3,ev4):
    return ev4 + exp(x) + 4*ev3 + x*ev4

def sol(y):

    h = 0.01
    rango = y - 1
    it = rango / 0.01
    ev0 = 2
    d = 1
    x1 = sp.diff(my_function, t)
    for i in range(int(it)):

        ev1 = x1.subs(x,ev0)
        ev2 = x2(d,ev0,ev1)
        ev3 = x3(d, ev1, ev2)
        ev4 = x4(d,ev2,ev3)
        ev5 = x5(d,ev3,ev4)
        ev0 = ev0 + h*ev1 + (h**2)/2 *ev2 + (h**3)/6 *ev3 + (h**4)/24 *ev4 +(h**5)/120 *ev5
        d += h
    return ev0

print(sol(3))
'''

t = sp.Symbol('t')
x = sp.Function('x')(t)


fun = x + sp.exp(t) + t*x

sols = taylor_method(5, fun, 2, 0.01, [1,3])
print(sols)
print(np.exp(9))
