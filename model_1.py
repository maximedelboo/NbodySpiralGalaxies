# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:42:49 2020

@author: Berend
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math       

def summation(i,n,m,x):
    out = 0
    r_i = x[i][0]
    
    for j in range(0,n):
        if i!=j:
            r_j = x[j][0]
            out = out + m[j]*(r_j - r_i)/(np.linalg.norm(r_j - r_i)**3)
    
    return out

# function that returns dx/dt
def model(x,t,n,G,m):    
    dxdt = np.array([[ [1,1] , [1,1] ]])
    x = np.reshape(x,(n,2,2))
    
    for i in range(0,n):
        v_i = x[i][1]
        
        dr_idt = v_i
        dv_idt = G*summation(i,n,m,x)
          
        dx_idt = np.array([[dr_idt, dv_idt]])
        dxdt = np.append(dxdt, dx_idt, axis = 0)
        
    dxdt = dxdt[1:]
    return np.reshape(dxdt, 4*n)

#number of bodies
n = 3

#masses sun, earth, mars kg
m =[1.989*(10**30), 5.972*(10**24), 6.39*(10**23)] 

#gravitational constant G
G = 6.67408*(10**(-11)) #SI

# initial condition
x0 = np.array([                         #r and v for all bodies
        [[0, 0],[0, 0]],                #sun 
        [[149600000000, 0],[0,30000]],  #earth
        [[218310000000, 0],[0,24000]]   #mars
        ])
    

t = np.linspace(0, 31556926, 3650)  #jaar in 10 stappen per dag

x0_1D = np.reshape(x0, 4*n)
  
sol = odeint(model, x0_1D,t, args = (n,G,m))
sol = np.reshape(sol,(3650,n,2,2))

plt.plot(sol[:,0,0,0],sol[:,0,0,1])
plt.plot(sol[:,1,0,0],sol[:,1,0,1])
plt.plot(sol[:,2,0,0],sol[:,2,0,1])

plt.show()


'''
# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()
'''