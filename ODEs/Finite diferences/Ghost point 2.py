import numpy as np
import matplotlib.pyplot as plt
from Ej7 import ghost_at_b

f = lambda x: -(np.pi**2) * np.cos(np.pi * x)
an = lambda x: np.cos(np.pi * x)

U, x = ghost_at_b(f,a=0,b=1/2,Ua=1,Uxb=-np.pi,n=50)

fig, ax = plt.subplots(1,1)
ax.plot(x, U, label='Solucion Numerica')
ax.plot(x, an(x), label='Solucion Analitica')
ax.grid()
ax.legend()
plt.show()