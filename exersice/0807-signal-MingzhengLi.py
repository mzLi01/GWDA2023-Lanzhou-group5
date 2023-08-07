# %%
import numpy as np
import matplotlib.pyplot as plt

a = [10, 3, 3]
A = 10

t = np.linspace(0, 1, 1000)
phi = sum([ai*t**(i+1) for i, ai in enumerate(a)])
s = A*np.sin(2*np.pi*phi)

plt.plot(t, s, marker='o', markersize=3)
plt.xlabel('t')
plt.ylabel('s')
