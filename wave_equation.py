

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
  Steps:
    1)  solve wave equation in 1d without source                    <==
    2)  include a source term
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

        self.x = []         # grid points
        self.h_p = []       # previous wave amplitude
        self.h = []         # current wave amplitude
        self.h_n = []       # next wave amplitude

        self.x = [ i*self.dx for i in range(N+1) ]
        self.h_p = [ self.f0(self.x[i]) for i in range(N+1) ]   # h(x,t=0) from initial waveform
        D = (self.c * self.dt / self.dx)**2
        for i in range(self.N+1) :                      # h(x,t=dt) from derivative of
            i_minus_1 = i - 1                           # of initial waveform
            i_plus_1  = i + 1

            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_1 = 0
            self.h.append(0.5*D*(self.h_p[i_plus_1] + self.h_p[i_minus_1]) + (1-D)*self.h_p[i] + self.dt*self.f0_prime(self.x[i]))

        self.h_n = [ 0 for i in range(N+1) ]

    
    def f0(self, x):                # initialize wave form
        self.x0 = self.L / 2.0      # starting position
        self.sigma = 0.05*self.L    # width
        k = math.pi / self.sigma
        gaussian = math.exp(-(x - self.x0)**2 / (2*self.sigma**2))
        return math.sin(k*(x-self.x0))*gaussian/10.


    def f0_prime(self, x):      # assumed the time dependence is only in sin factor
        self.x0 = self.L / 2.0
        self.sigma = 0.05*self.L
        k = math.pi / self.sigma
        gaussian = math.exp(-(x - self.x0)**2 / (2*self.sigma**2))
        return math.cos(k*(x - self.x0))*gaussian*self.c / (10.*k)


    def FTCS(self):             # not really FTCS method but I didn't have another name
        for i in range(self.N+1):
            i_minus_1 = i - 1
            i_plus_1 = i + 1

            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_1  = 0
            
            D = (self.c * self.dt / self.dx )**2
            self.h_n[i] = -self.h_p[i] + 2.*self.h[i] + D*(self.h[i_plus_1] + self.h[i_minus_1] - 2.*self.h[i])


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
