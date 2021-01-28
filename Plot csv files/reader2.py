import numpy as np
import sys
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers

from matplotlib import rcParams

rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'  # ffmpeg path voor videos

fig, ax = plt.subplots(figsize=(8, 8))

filename = 'C:/Users/tip/Files/Git projects/Planeten-enal/huts.csv'
# path van de csv, best bij de python file
dat = np.genfromtxt(filename, delimiter=';')
dat = dat[:, :-1]
row1 = dat[0]
poii = np.random.randint(int(len(row1)/2))
poii2 = np.random.randint(int(len(row1)/2))
poii3 = np.random.randint(int(len(row1)/2))
poii4 = np.random.randint(int(len(row1)/2))
row1 = row1.reshape(int(len(row1) / 2), 2).transpose()
rs = 1.2
ax.set(xlim=(min(row1[0]) * rs, max(row1[0]) * rs), ylim=(min(row1[1]) * rs, max(row1[1]) * rs))

line = ax.plot(row1[0], row1[1], '.', ms=1)[0]

ms_line = 10

line2 = ax.plot(row1[0][poii], row1[1][poii], '.', ms=ms_line, c='red')[0]
line3 = ax.plot(row1[0][poii2], row1[1][poii2], '.', ms=ms_line, c='red')[0]
line4 = ax.plot(row1[0][poii3], row1[1][poii3], '.', ms=ms_line, c='red')[0]
line5 = ax.plot(row1[0][poii4], row1[1][poii4], '.', ms=ms_line, c='red')[0]
#line4 = ax.plot(row1[0][900], row1[1][900], '.', ms=10, c='red')[0]


def animate(i):
    ax.set_title(str(i))
    row = dat[(i % len(dat))]
    row = row.reshape(int(len(row) / 2), 2).transpose()
    line.set_ydata(row[1])
    line.set_xdata(row[0])
    line2.set_ydata(row[1][poii])
    line2.set_xdata(row[0][poii])
    line3.set_ydata(row[1][poii2])
    line3.set_xdata(row[0][poii2])
    line4.set_ydata(row[1][poii3])
    line4.set_xdata(row[0][poii3])
    line5.set_ydata(row[1][poii4])
    line5.set_xdata(row[0][poii4])


anim = FuncAnimation(fig, animate, interval=10,
                     frames=5000)  # je kan rondkloten met frames en interval voor de snelheid
plt.draw()
plt.title(filename)
plt.show()

# VOOR VIDEO

# from matplotlib import rcParams
# rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\bin\ffmpeg.exe' #ffmpeg path voor videos

# Writer = writers['ffmpeg']
# writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

# anim.save('barnes-HUTSpar.mp4', writer=writer)
