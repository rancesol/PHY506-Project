

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
    def __init__(self, method='Lax_Wendroff()', N=100, L=100, c=1.0):
        self.L = L                  # Length of system
        self.N = N                  # number of cells
        self.dx = float(L)/float(N) # cell size
        self.c = c                  # speed of waves
        self.t = 0.0                # time
        self.dt = 0.1*self.dx #/c         # time step size
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function

        self.x = []         # grid points
        self.h0 = []        # initial wave form
#        self.h_p = []       # previous wave amplitude
        self.h = []         # current wave amplitude
        self.h_n = []       # next wave amplitude
        self.h_nn = []      # wave amp. at next-next time step (due to 2nd order derivatives)

        self.dx = L / float(N)
        self.x = [ i*self.dx for i in range(N+1) ]
        self.h0 = [ self.f0(self.x[i]) for i in range(N+1) ]
#        self.h_p = [ self.h0[i] for i in range(N+1) ]
        self.h = [ self.h0[i] for i in range(N+1) ]
        self.h_n = [ self.h[i] for i in range(N+1) ]
        self.h_nn = [ self.h_n[i] for i in range(N+1) ]

    
    def f0(self, x):            # initialize wave form
        self.x0 = self.L / 2.0  # starting position
        self.sigma = 0.1*self.L # width

        k = math.pi / self.sigma
        gaussian = math.exp(-(x - self.x0)**2 / (2*self.sigma**2))
        return math.cos(k*(x-self.x0))*gaussian


    def FTCS(self):
        for i in range(self.N+1):
            i_minus_1 = i - 1
            i_plus_1 = i + 1
            i_plus_2 = i + 2

            if i == self.N - 1 :
                i_plus_2 = 0
            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_2 = 1
                i_plus_1 = 0
            
            self.h_nn[i] = 2*self.h_n[i] - self.h[i] + ( self.h[i_plus_2] - 2*self.h[i_plus_1] + self.h[i] ) * self.c**2 * self.dt**2 / self.dx**2
#            D = (self.c * self.dt / self.dx)**2
#            self.h_n[i] = self.h_p[i] - 2*self.h[i] + D*(self.h[i_plus_1] + self.h[i_minus_1] - 2*self.h[i])
            #self.h_nn[0] = self.h_nn[self.N]

    def Lax_Wendroff(self) :
        D = (self.c * self.dt / self.dx)**2 / 2.0
        for i in range(self.N+1):
            i_plus_1  = i + 1
            i_plus_2  = i + 2
            i_minus_1 = i - 1

            if i == self.N - 1 :
                i_plus_2 = 0
            if i == self.N :
                i_plus_2 = 1
                i_plus_1 = 0
            if i == 0 :
                i_minus_1 = self.N
            
            self.h_nn[i] = 2*self.h_n[i] + self.h[i] + ( self.h[i_plus_2] - 2*self.h[i_plus_1] - self.h[i] ) * self.dt**2 / self.dx**2
            self.h_nn[i] += D * (self.h[i_plus_1] + self.h[i_minus_1] - 2 * self.h[i])


    def take_step(self):
        eval('self.' + self.method )
        swap = self.h
        self.h = self.h_n
        self.h_n = self.h_nn
        self.h_nn = swap

#        swap = self.h
#        self.h = self.h_n
#        self.h_n = self.h_p
#        self.h_p = swap

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
        self.ax.set_ylim(-2.,2.)
        initvals = [ ix for ix in self.dalembertian.h]
        self.line, = self.ax.plot(initvals)


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


dalembertian = dAlembertian( method='FTCS()', N=100, L=100, c=1. )
animator = Animator( dalembertian=dalembertian )
animator.animate()
plt.show()
