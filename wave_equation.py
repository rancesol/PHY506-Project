

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
  Steps:
    1)  solve wave equation in 1d without source
    2)  include a source term                                       <==
    3)  make source term an infall under the slow motion condition
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

class dAlembertian :
    def __init__(self, method='FTCS', N=100, L=100, c=1.0):
        self.L = L                  # Length of system
        self.N = N                  # number of cells
        self.dx = float(L)/float(N) # cell size
        self.c = c                  # speed of waves
        self.G = 6.67e-11           # Gravitational constant
        self.t = 0.0                # time
        self.dt = 0.001              # time step size
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function

        self.Tloc = self.L / 20.0   # source location
        self.omega_o = 1.0          # initial orbital frequency of source masses
        self.w = []                 # orbital frequency of source masses
        self.M = 1.0                # mass of source objects (taken to be the same)
        self.R_o = 1.0              # radius of orbit
        self.R_M = 1.0               # radius of masses

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
        Tmax  = 500.            # strength of source
        if abs(x-self.Tloc) <= self.dx :
            self.omega()        # update frequency
            if abs(self.w[-1] - math.sqrt(self.G*self.M/(4.*self.R_M**3))) < 0.1 :
                    return 0
            else :
                    return Tmax * math.sin(2.*self.w[-1]*t)
        else :
            return 0


    def omega(self) :
        E_i = - self.G*self.M**2 / (4.*self.R_o)
        L   = (8./5.)*((2.*self.G)**(4./3.))*(self.M**(10./3.))*(self.w[-1]**(10./3.))
        self.w.append(self.w[-1] - (3./2.)*self.w[-1]*L*self.dt/(E_i-L))
        omega_file.write(repr(self.t) + '\t' + repr(2.*self.w[-1]) + '\n')


    def FTCS(self):             # not really FTCS method but I didn't have another name
        for i in range(self.N+1):
            i_minus_1 = i - 1
            i_plus_1 = i + 1

            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_1  = 0
            
            D = (self.c * self.dt / self.dx )**2
            self.h_n[i] = -self.h_p[i] + 2.*self.h[i] + D*(self.h[i_plus_1] + self.h[i_minus_1] - 2.*self.h[i]) + self.T(self.dx*i, self.t) * self.dt**2

            if self.x[i] < self.Tloc :
                self.h_n[i] *= (math.exp(math.log(2)*self.x[i]/self.Tloc) - 1)


    def take_step(self):
        eval('self.' + self.method )
        swap = self.h
        self.h = self.h_n
        self.h_n = self.h_p
        self.h_p = swap

        self.t += self.dt
        self.step_number += 1

    
    def sav_h(self, plot_number):
        file_name = 'h_' + repr(plot_number) + '.data'
        file = open(file_name, 'w')
        for i in range(self.N):
            file.write(repr(self.x[i]) + '\t' + repr(self.h[i]) + '\n')
        file.close()
        print ' saved h(x,t) at t = ', self.t, ' in ', file_name


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


dalembertian = dAlembertian( method='FTCS()', N=10000, L=1000, c=50. )
animator = Animator( dalembertian=dalembertian )
animator.animate()
plt.show()
omega_file.close()
