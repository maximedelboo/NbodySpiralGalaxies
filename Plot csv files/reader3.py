import numpy as np
import sys
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers

from matplotlib import rcParams

rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'  # ffmpeg path voor videos

fig, ax = plt.subplots(figsize=(8, 8))

filename = 'C:/Users/tip/Files/Git projects/Planeten-enal/huts b.csv'
# path van de csv, best bij de python file
dat = np.genfromtxt(filename, delimiter=';')
dat = dat[:, :-1]
row1 = dat[0]
poii = np.random.randint(int(len(row1)/2))
poii2 = np.random.randint(int(len(row1)/2))
#poii = 600
#poii2 = 8000
#poii3 = 4400
#poii4 = 4100
poii3 = np.random.randint(int(len(row1)/2))
poii4 = np.random.randint(int(len(row1)/2))
row1 = row1.reshape(int(len(row1) / 2), 2).transpose()
rs = 1.2
ax.set(xlim=(min(row1[0]) * rs, max(row1[0]) * rs), ylim=(min(row1[1]) * rs, max(row1[1]) * rs))

line = ax.plot(row1[0], row1[1], '.', ms=1)[0]

ms_line = 10

print(poii)
print(poii2)
print(poii3)
print(poii4)

Gvar = 6.67408e-11
m_sol = 1.98847e30
m_center_mass = 700
radius_gal = 1e18
omega = np.sqrt(Gvar*m_center_mass*m_sol / (radius_gal / 2)**3)
dt = 35e12

line2 = ax.plot(row1[0][poii], row1[1][poii], ms=10, c='red')[0]
line3 = ax.plot(row1[0][poii2], row1[1][poii2], ms=10, c='green')[0]
line4 = ax.plot(row1[0][poii3], row1[1][poii3], ms=10, c='yellow')[0]
line5 = ax.plot(row1[0][poii4], row1[1][poii4], ms=10, c='black')[0]
line6 = ax.plot([0, radius_gal], [0, 0], ms=5, c='black')[0]



p_c = np.linspace(0, 1, 100) * 2 * np.pi
r_in = (Gvar * m_center_mass * m_sol/(4*omega**2))**(1/3)
r_out = (9*Gvar * m_center_mass * m_sol/(4*omega**2))**(1/3)
x_c_i = np.cos(p_c)*r_in
y_c_i = np.sin(p_c)*r_in
x_c_o = np.cos(p_c) * r_out
y_c_o = np.sin(p_c) * r_out

line7 = ax.plot(x_c_i, y_c_i, ms=5, c='purple')[0]
line8 = ax.plot(x_c_o, y_c_o, ms=5, c='purple')[0]
#line4 = ax.plot(row1[0][900], row1[1][900], '.', ms=10, c='red')[0]

p1x = []
p2x = []
p3x = []
p4x = []
p1y = []
p2y = []
p3y = []
p4y = []


def animate(i):
    t = i*dt
    T = omega*t
    #x_a = np.linspace(0, radius_gal, 100)
    #x_y = np.linspace(0, radius_gal, 100)
    o_x = radius_gal * np.cos(T)
    o_y = radius_gal * np.sin(T)

    ax.set_title(str(i))
    row = dat[(i % len(dat))]
    row = row.reshape(int(len(row) / 2), 2).transpose()
    p1y.append(row[1][poii])
    p1x.append(row[0][poii])
    p2y.append(row[1][poii2])
    p2x.append(row[0][poii2])
    p3y.append(row[1][poii3])
    p3x.append(row[0][poii3])
    p4y.append(row[1][poii4])
    p4x.append(row[0][poii4])
    s_i = min(len(p4x), 200)
    line.set_ydata(row[1])
    line.set_xdata(row[0])
    line2.set_ydata(p1y[-s_i:])
    line2.set_xdata(p1x[-s_i:])
    line3.set_ydata(p2y[-s_i:])
    line3.set_xdata(p2x[-s_i:])
    line4.set_ydata(p3y[-s_i:])
    line4.set_xdata(p3x[-s_i:])
    line5.set_ydata(p4y[-s_i:])
    line5.set_xdata(p4x[-s_i:])

    line6.set_xdata([0, o_x])
    line6.set_ydata([0, o_y])


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
