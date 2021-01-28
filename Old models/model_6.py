# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 13:11:57 2020

@author: Berend
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import pickle

"""
Create Your Own N-body Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz
Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity
"""
@jit
def getAcc( pos, mass, G, softening, dt, i ):
    """
    Calculate the acceleration on each particle due to Newton's Law 
    pos  is an N x 2 matrix of positions
    mass is an N x 1 vector of masses
    G is Newton's Gravitational constant
    softening is the softening length
    a is N x 2 matrix of accelerations
    """

    # positions r = [x,y,z] for all particles
    x = pos[:,0:1]
    y = pos[:,1:2]
    
    # matrix that stores all pairwise particle separations: r_j - r_i
    dx = x.T - x
    dy = y.T - y
    
    # matrix that stores 1/r^3 for all particle pairwise particle separations 
    inv_r3 = (dx**2 + dy**2 + softening**2)
    inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)
    
    #pertubator
    cst = 1e23
    
    n = (x**2 +y**2 +softening**2)
    
    omega = 2*np.pi/(150*dt)
    t = i*dt
    
    T = omega*t
    
    AX = cst*( x*(x**2-y**2)*np.cos(T) + 2*y*(x**2)*np.sin(T) - 4*(y**2)*x*np.cos(T) + 2*y*(x**2-y**2)*np.sin(T) )*(n**-2.5)
    AY = cst*( y*(x**2-y**2)*np.cos(T) + 2*x*(y**2)*np.sin(T) + 4*(x**2)*y*np.cos(T) - 2*x*(x**2-y**2)*np.sin(T) )*(n**-2.5)
    
    print('T is:')
    print(T)
    
    '''
    r^2 cos(2 theta):
    AX = cst*( (4*x*y**2) * n**(-2.5) + (x*y**2 - x**3)*n**(-1.5) ) 
    AY = cst*( (y**3 - y*x**2)*n**(-1.5) - (4*y*x**2)*n**(-2.5))
    
    (r)^-1 cos(2 theta):
    AX = cst*( x*(x**2 - 5*y**2) * n**(-2.5))
    AY = cst*( y*(5*x**2 - y**2) * n**(-2.5))
    
    '''
    
    #central mass
    inv_r3_M = (x**2 + y**2 + softening**2)**(-1.5)
    M =  1e4*1.989*(10**30)
    
    
    #accelerations
    
    ax = G * (dx * inv_r3) @ mass - G * (x * inv_r3_M) * M + AX
    ay = G * (dy * inv_r3) @ mass - G * (y * inv_r3_M) * M + AY
          	
    # pack together the acceleration components
    a = np.hstack((ax,ay))
    
    return a
@jit
def main():
    """ N-body simulation """
    	
    # Simulation parameters
    N         = 2000    # Number of particles
    t         = 0      # current time of the simulation
    tEnd      = 5*3.5e7   # time at which simulation ends
    dt        = 0.5*3.5e4   # timestep
    rad_gal = 1e12     #radius of the galaxy
    softening = 0.1* rad_gal   # softening length
    G         = 6.67e-11    # Newton's Gravitational Constant
    plotRealTime = True # switch on for plotting as the simulation goes along
    
    	
    # Generate Initial Conditions
    np.random.seed(17)            # set the random number generator seed
    
    #masses
    m_sol = 1.989*(10**30)      #solar mass    	
    mass = m_sol * np.ones((N,1))  # N stars with mass of the sun
    
    #positions
    phi = 2*np.pi*np.random.rand(N)     #phi is the array with azimuthal angles 
    r = np.arange(501,N+501)*((rad_gal)/N)    #r is the array with the radii
    pos  =  np.transpose(np.array([r*np.cos(phi), r*np.sin(phi)]))   #positions from r and phi
    pos = pickle.load( open("C:\\Users\\Berend\\pos.pickle", "rb" ) )
    
    #velocities
    v = np.sqrt( G*m_sol*1e4/r)  #the magnitudes of the velocity vectors based on r
    vel  = np.transpose([v*(-np.sin(phi)) , v*np.cos(phi)])
    vel  = pickle.load( open("C:\\Users\\Berend\\vel.pickle", "rb" ) ) 
    
    	
    # Convert to Center-of-Mass frame
    #vel -= np.mean(mass * vel,0) / np.mean(mass)

    	
    # calculate initial gravitational accelerations
    acc = getAcc( pos, mass, G, softening, dt, 0 )
  	
    # number of timesteps
    Nt = int(np.ceil(tEnd/dt))
    	
    # save energies, particle orbits for plotting trails
    pos_save = np.zeros((N,2,Nt+1))
    pos_save[:,:,0] = pos
    t_all = np.arange(Nt+1)*dt
    	
    # prep figure
    fig = plt.figure(figsize=(4,5), dpi=80)
    grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
    ax1 = plt.subplot(grid[0:2,0])
    
    # Simulation Main Loop
    for i in range(Nt):
        # (1/2) kick
        vel += acc * dt/2.0
        	
        # drift
        pos += vel * dt
        		
        # update accelerations
        acc = getAcc( pos, mass, G, softening, dt, i )
        		
        # (1/2) kick
        vel += acc * dt/2.0
        
        # Convert to Center-of-Mass frame
        #pos -= np.mean(mass * pos,0) / np.mean(mass)
        		
        # update time
        t += dt
        		
        # save energies, positions for plotting trail
        pos_save[:,:,i+1] = pos
        
        print('i is:')
        
        print(i)
        '''
        if i == 500:
            pickle.dump( pos, open( "C:\\Users\\Berend\\pos.pickle", "wb" ))
            pickle.dump( vel, open( "C:\\Users\\Berend\\vel.pickle", "wb" ))
        '''
        if i == 1000:
            pickle.dump( pos, open( "C:\\Users\\Berend\\pos_hm.pickle", "wb" ))
        	
        # plot in real time
        if plotRealTime or (i == Nt-1):
            plt.sca(ax1)
            plt.cla()
            xx = pos_save[:,0,max(i-50,0):i+1]
            yy = pos_save[:,1,max(i-50,0):i+1]
            #plt.scatter(xx,yy,s=1,color=[.7,.7,1])          #draws piece of the path
            plt.scatter(pos[:,0],pos[:,1],s=10,color='blue')
            ax1.set(xlim=(-2, 2), ylim=(-2, 2))
            ax1.set_aspect('equal', 'box')
            
            
            ax1.set_xticks([-3.5*rad_gal, -0.5*rad_gal, 0, 0.5*rad_gal, 3.5*rad_gal])
            ax1.set_yticks([-3.5*rad_gal, -0.5*rad_gal, 0, 0.5*rad_gal, 3.5*rad_gal])
            		
            plt.pause(0.001)
            
        
    	    
    	
    	
    # add labels/legend
    plt.xlabel('x')
    plt.ylabel('y')
    
    # Save figure
    plt.savefig('nbody.png',dpi=240)
    plt.show()
    	    
    return 0
	  
if __name__== "__main__":
    main()