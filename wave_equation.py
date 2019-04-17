

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


class dAlembertian :
    def __init__(self, method='FTCS', N=100, L=100, c=1.0):
        self.L = L                  # Length of system
        self.N = N                  # number of cells
        self.dx = float(L)/float(N) # cell size
        self.c = c                  # speed of waves
        self.t = 0.0                # time
        self.dt = 0.01              # time step size
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function

        self.Tloc = self.L / 20.0    # source location

        self.x = []         # grid points
        self.h_p = []       # previous wave amplitude
        self.h = []         # current wave amplitude
        self.h_n = []       # next wave amplitude

        self.x = [ i*self.dx for i in range(N+1) ]
        self.h_p = [ 0 for i in range(N+1) ]    # take initial waveform as zero everywhere
        self.h   = [ 0 for i in range(N+1) ]
        self.h_n = [ 0 for i in range(N+1) ]


    def T(self, x, t):          # source term
        omega = 1.*(2.*math.pi) # oscillation frequency
        Tmax  = 30.             # strength of source
        if abs(x-self.Tloc) <= self.dx :
            return Tmax * math.cos(omega*t**2)
        else :
            return 0


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


dalembertian = dAlembertian( method='FTCS()', N=500, L=100, c=10. )
animator = Animator( dalembertian=dalembertian )
animator.animate()
plt.show()
