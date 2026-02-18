import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd


try:
    data = np.loadtxt('poukladane_dane.csv')
except OSError:
    print("Nie znaleziono pliku 'wyniki.txt'")


# print(data[0])

fig, ax = plt.subplots(figsize=(6, 4))
y1 = data[1]
x = np.linspace(0.4, 2.0, 32)
y2 = data[2]
y3 = data[3]
y4 = data[4]
y5 = data[5]
y6 = data[6]
y = data[0]

ax.plot(x, y)
ax.plot(x, y1)
ax.plot(x, y2)
ax.plot(x, y3)
ax.plot(x, y4)
ax.plot(x, y5)
ax.plot(x, y6)
ax.set_title("Energia")
ax.set_xlabel("alpha")
ax.set_ylabel("E")

plt.savefig("Energia.png")