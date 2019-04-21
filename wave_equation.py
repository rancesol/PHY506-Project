

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
  Steps:
    1)  solve wave equation in 1d without source
    2)  include a source term
    3)  make source term an infall under the slow motion condition <==
    4)  perform spectral analysis and get the "chirp"
    5)  remove slow motion condition
"""
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


import math
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

omega_file = open("omega.data", "w")

# Geometrized units:
# G = 6.67e-11 => 1
# c = 3e8 => 1
# 1M_sun = 2e30kg => 1.48e3 m
# 1s = 3e8 m
# 1pc = 3e16m



class dAlembertian :
    def __init__(self, method='KDRS', N=100, L=100, c=1.0):
        self.L = L                  # Length of system
        self.N = N                  # number of cells
        self.dx = float(L)/float(N) # cell size
        self.c = c                  # speed of waves
        self.G = 6.67e-11           # Gravitational constant
        self.Msol = 2.e30           # solar mass in meters
        self.t = 0.0                # time
        self.dt = 0.01              # time step size
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function

        self.Tloc = self.L / 50.0   # source location
        self.w = []                 # orbital frequency of source masses
        self.M = 32.5*self.Msol     # mass of source objects (taken to be the same)
        self.R_o = 1.e11            # starting orbital radius (about 1AU)
        self.R_M = 350.e3           # radius of masses (350km -- approx radius of BHs in GW150914)
        self.omega_o = math.sqrt(self.G*self.M/(4.*self.R_o**3))    # initial freq.
        self.omega_max = math.sqrt(self.G*self.M/(4.*self.R_M**3))  # freq. at merger

        self.x = []         # grid points
        self.h_p = []       # previous wave amplitude
        self.h = []         # current wave amplitude
        self.h_n = []       # next wave amplitude

        # initialize
        self.w.append(self.omega_o)
        self.x   = [ i*self.dx for i in range(N+1) ]
        self.h_p = [ 0 for i in range(N+1) ]    # take initial waveform as zero everywhere
        self.h   = [ 0 for i in range(N+1) ]
        self.h_n = [ 0 for i in range(N+1) ]


    def T(self, x, t):          # source term
        Tmax  = 500.            # strength of source (kept const. for simplicity)
        if abs(x-self.Tloc) < self.dx/2. :
            self.w_update()     # update frequency
            if self.w[-1] >= self.omega_max :
                    return 0
            else :
                    return Tmax * math.sin(2.*self.w[-1]*t) # GWs have 2*freq. of source
        else :
            return 0


    def w_update(self) :    # update orbital frequency
        E_i = - self.G*self.M**2 / (4.*self.R_o)
        E_GW = 0.           # energy radiated away in GWs
        E_GW += (8./5.)*((2.*self.G)**(4./3.))*(self.M**(10./3.))*(self.w[-1]**(10./3.))*self.dt
        self.w.append(self.w[-1] - (3./2.)*self.w[-1]*E_GW*self.dt/(E_i-E_GW))
        if self.w[-1] >= self.omega_max :
            omega_file.write(repr(self.t) + '\t' + repr(0.0) + '\n')
        else :
            omega_file.write(repr(self.t) + '\t' + repr(2.*self.w[-1]) + '\n')


    def KDRS(self):             # KDRS method (not a thing!)
        for i in range(self.N+1):
            i_minus_1 = i - 1
            i_plus_1 = i + 1

            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_1  = 0
            
            D = (self.c * self.dt / self.dx )**2
            self.h_n[i] = -self.h_p[i] + 2.*self.h[i] + D*(self.h[i_plus_1] +
                    self.h[i_minus_1] - 2.*self.h[i]) + self.T(self.dx*i, self.t) * self.dt**2

            if self.x[i] < self.Tloc :  # kill off everything leftward of the source
                self.h_n[i] *= 0.0


    def take_step(self):
        eval('self.' + self.method )
        swap = self.h
        self.h = self.h_n
        self.h_n = self.h_p
        self.h_p = swap

        self.t += self.dt
        self.step_number += 1



class Animator :
    def __init__(self, periodic = True, dalembertian = None):
        self.avg_times = []
        self.dalembertian = dalembertian
        self.t = 0.
        self.fig, self.ax = plt.subplots()
        self.ax.set_ylim(-.5,.5)
        initvals = [ ix for ix in self.dalembertian.h]
        x = [ ix for ix in self.dalembertian.x ]
        self.line, = self.ax.plot(x,initvals)


    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,

    
    def time_step(self) :
        dalembertian.take_step()
        yield [ix for ix in self.dalembertian.h]


    def create_widgets(self) :
        self.QUIT = Button(self, text='QUIT', command=self.quit)
        self.QUIT.pack(side=BOTTOM)

        self.draw = Canbas(self, width='600', height='400')
        self.draw.pack(side=TOP)


    def animate(self) :
        self.ani = animation.FuncAnimation( self.fig,
                                            self.update,
                                            self.time_step,
                                            interval=50,
                                            blit=False )


dalembertian = dAlembertian( method='KDRS()', N=10000, L=1000, c=10. )
animator = Animator( dalembertian=dalembertian )
animator.animate()
plt.show()
omega_file.close()
