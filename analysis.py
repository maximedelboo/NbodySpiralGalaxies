# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 15:39:59 2020

@author: Berend
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

positions = pickle.load( open("C:\\Users\\Berend\\pos_hm.pickle", "rb" ) )
velocities = pickle.load( open("C:\\Users\\Berend\\vel.pickle", "rb" ))


x = np.linspace(-0.5e12, 0.5e12, 100)
y = np.linspace(-0.5e12, 0.5e12, 100)

x,y = np.meshgrid(x,y)

def force(x,y):
    cst = 1e23
    
    n = (x**2 +y**2)
    
    T = np.pi/4
    
    AX = cst*( x*(x**2-y**2)*np.cos(T) + 2*y*(x**2)*np.sin(T) - 4*(y**2)*x*np.cos(T) + 2*y*(x**2-y**2)*np.sin(T) )*(n**-2.5)
    AY = cst*( y*(x**2-y**2)*np.cos(T) + 2*x*(y**2)*np.sin(T) + 4*(x**2)*y*np.cos(T) - 2*x*(x**2-y**2)*np.sin(T) )*(n**-2.5)
    
    return AX**2 + AY**2

def density(x,y,pos):
    epsilon = 5e10
    dx = np.abs(x-pos[:,0])
    dy = np.abs(y-pos[:,1])
    
    d = np.sqrt(dx**2 + dy**2)
    
    return len(d[ d < epsilon ])

def window(a,lwbd, upbd):
    if a<upbd and a>lwbd:
        return a
    elif a >= upbd:
        return upbd
    
    else:
        return lwbd
        

density_vec = np.vectorize(density, excluded = ['pos'])
window_vec = np.vectorize(window, excluded = ['lwbd', 'upbd'])
#force_vec = np.vectorize(force)

z = density_vec(x,y,pos = positions)
#z = force_vec(x,y)
z = window_vec(z, lwbd = 0, upbd = 100)

ax = sb.heatmap(z)
ax.invert_yaxis()
fig =ax.get_figure()
fig.savefig('heatmap.png')

plt.plot(np.sqrt(positions[:,0]**2 + positions[:,1]**2),np.sqrt(velocities[:,0]**2 + velocities[:,1]**2), 'o')