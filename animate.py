import matplotlib.animation as animation
from matplotlib import style
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

def animate(i):
    r_x = np.loadtxt('C:/Users/tip/Files/Git projects/modeling2b/model/model_x.txt')
    r_y = np.loadtxt('C:/Users/tip/Files/Git projects/modeling2b/model/model_y.txt')
    ax.clear()
    ax.plot(r_x, r_y, ls='', marker='.')
    ax.set_aspect('equal')


ani = animation.FuncAnimation(fig, animate, interval=500)

#ax.plot(r_x, r_y, ls='', marker='.')

#ax.set_aspect('equal')

plt.show()

