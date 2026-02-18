import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd


try:
    data = np.loadtxt('wekt_wlas.dat')
except OSError:
    print("Nie znaleziono pliku 'wyniki.txt'")


print(data[:,0])

fig, ax = plt.subplots(figsize=(6, 4))
y1 = data[:,1]
x = data[:,0]
y2 = data[:,2]
y3 = data[:,3]
y4 = data[:,4]
# y5 = data[:,5]
# y6 = data[:,6]

ax.plot(x, y1)
ax.plot(x, y2)
ax.plot(x, y3)
ax.plot(x, y4)
# ax.plot(x, y5)
# ax.plot(x, y6)
ax.set_title("Psi")
ax.set_xlabel("x")
ax.set_ylabel("Psi")

plt.savefig("psi.png")