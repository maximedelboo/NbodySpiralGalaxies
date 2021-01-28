import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,writers


from matplotlib import rcParams
rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\bin\ffmpeg.exe' #ffmpeg path voor videos


fig, ax = plt.subplots(figsize=(8,8))

filename = 'huts.csv' #path van de csv, best bij de python file
dat = np.genfromtxt(filename, delimiter=';')
dat = dat[:,:-1]
row1 = dat[0]
row1 = row1.reshape(int(len(row1)/2),2).transpose()
rs = 1.2
ax.set(xlim = (min(row1[0])*rs,max(row1[0])*rs),ylim = (min(row1[1])*rs,max(row1[1])*rs))

line  = ax.plot(row1[0],row1[1], '.', ms =1)[0]

def animate(i):
    row = dat[(i % len(dat))]
    row = row.reshape(int(len(row)/2),2).transpose()
    line.set_ydata(row[1])
    line.set_xdata(row[0])
    

anim = FuncAnimation(fig,animate,interval=25,frames = 1000) #je kan rondkloten met frames en interval voor de snelheid
plt.draw()
plt.title(filename)
plt.show()


#VOOR VIDEO

#from matplotlib import rcParams
#rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\bin\ffmpeg.exe' #ffmpeg path voor videos

#Writer = writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

#anim.save('barnes-HUTSpar.mp4', writer=writer)

