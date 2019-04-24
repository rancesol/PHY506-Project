

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
    def __init__(self, method='KDRS', N=100, L=100):
        self.L = L                  # Length of system
        self.N = N                  # number of cells
        self.dx = float(L)/float(N) # cell size
        self.c = 2.99792e8 #c                  # speed of waves
        self.G = 6.67408e-11           # Gravitational constant
        self.Msol = 1.9891e30           # solar mass in meters
        self.t = 0.0                # time
        self.dt = 0.01              # time step size
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function

        self.Tloc = self.L / 50.0   # source location
        self.w = []                 # orbital frequency of source masses
        self.M = 32.5*self.Msol     # mass of source objects (taken to be the same)
        self.eta = 1./4             # symmetric mass difference
        self.R_o = 2.e6            # starting orbital radius (about 1AU)
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
            if self.w[-1] >= self.omega_max :
                return 0
            else :
                self.w_update4()    # update freq.
                print self.w[-1]
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


    def w_update2(self) :   # update orbit freq. (Eq.2.33 Analysis of GW Data - Jaranowski, Krolak)
        self.w.append(self.omega_o*(1. -
            256.*(self.G*self.M)**(5./3)*self.omega_o**(8./3)*self.t/(5.*self.c**5))**(-3./8))
        if self.w[-1] >= self.omega_max :
            omega_file.write(repr(self.t) + '\t' + repr(0.0) + '\n')
        else :
            omega_file.write(repr(self.t) + '\t' + repr(2.*self.w[-1]) + '\n')

    
    def w_update3(self) :
        eta = self.eta
        x_o = (self.G*self.M*self.omega_o)**(2./3)/self.c**2
        t_hat = self.c**3*self.t/(self.G*self.M)
        tau = eta*(self.F(x_o) - t_hat)/5
        x = self.Y(tau)
        self.w.append(x**(3./2)*self.c**3/(self.G*self.M))

        if self.w[-1] >= self.omega_max :
            omega_file.write(repr(self.t) + '\t' + repr(0.0) + '\n')
        else :
            omega_file.write(repr(self.t) + '\t' + repr(2.*self.w[-1]) + '\n')


    def w_update4(self) :
        t_c = 5.*self.c**5*self.R_o**4/(256.*self.G**3*2*self.M**3)
        self.w.append(2.**(1./8)*5.**(3./8)*self.c**(15./8)*(t_c - self.t)**(-3./8)/(8.*math.pi*self.G**(5./8)*self.M**(5./8)))
        if self.w[-1] >= self.omega_max :
            omega_file.write(repr(self.t) + '\t' + repr(0.0) + '\n')
        else :
            omega_file.write(repr(self.t) + '\t' + repr(2.*self.w[-1]) + '\n')


    def F(self, x) :
        eta = self.eta
        gamma_E = 0.577216
        F  = 1. + (743./252 + eta*11./3)*x - (32./5)*math.pi*x**(3./2)
        F += (3058673./508032 + eta*5429./504 + eta**2*617./72)*x**2
        F += (-7729./252 + eta*13./3)*math.pi*x**(5./2)
        F += (110052469856691./23471078400 + 128.*math.pi**2/3 + 6848.*gamma_E/105 + 3424.*math.log(16.*x)/105)*x**3
        F += ((3147553127./3048192 - 451.*math.pi**2/12)*eta - 15211.*eta**2/1728 + 25565*eta**3/1296)*x**3
        F += (-15419335./127008 - 75703.*eta/756 + 14809*eta**2/378)*math.pi*x**(7./2)
        F *= 5./(256.*eta*x**4)
        return F

    def Y(self, tau) :
        eta = self.eta
        gamma_E = 0.577216
        Y  = 1. + (743./4032 + 11.*eta/48)*tau**(-1./4) - math.pi*tau**(-3./8)/5
        Y += (19583./254016 + 24401.*eta/193536 + 31.*eta**2/288)*tau**(-1./2)
        Y += (-11891./53760 + 109.*eta/1920)*math.pi*tau**(-5./8)
        Y += (math.pi**2/6 - 10052469856691./6008596070400 + 107.*gamma_E/420 - 107.*math.log(tau/256)/3360)*tau**(-3./4)
        Y += (-113868647./433520640 - 31821.*eta/143360 + 294941.*eta**2/3870720)*math.pi*tau**(-7./8)
        Y *= tau**(-1./4)/4
        return Y


    def KDRS(self):             # KDRS method (not a thing!)
        for i in range(self.N+1):
            i_minus_1 = i - 1
            i_plus_1 = i + 1

            if i == 0 :
                i_minus_1 = self.N
            if i == self.N :
                i_plus_1  = 0
            
            D = 1. #(self.c * self.dt / self.dx )**2
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
                                            interval=1,
                                            blit=False )


dalembertian = dAlembertian( method='KDRS()', N=10000, L=1000 )
animator = Animator( dalembertian=dalembertian )
animator.animate()
plt.show()
omega_file.close()
