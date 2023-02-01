import scipy.constants as sc
import numpy as np
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  Parameters  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c = sc.c # speed of light

# courant_factor condition < 1
courant_factor = 0.4
#System 1D
N_c = 400 # cells number
dx = sc.parsec * 6.25
L = N_c*dx # total distance


# Time 
N_t = 800 #iteration number
dt = dx*courant_factor / c
tf = N_t*dt # final time

# Lists for Numerical Integration
N = 1e-20*np.ones(N_c)
F = np.zeros(N_c)
F_GLF = np.zeros(N_c)
P = np.zeros(N_c)
P_GLF = np.zeros(N_c)
Khi = np.zeros(N_c)

# List for Ionization
x = np.ones(N_c)*0.0012

# Physical parameters
T = 2e4 # Temperature in Kelvin
rho = 1e-2 # Number density of hydrogen nuclei in cm^-3, must stay in cm
sigma = 1.63*10**-18 #cm^-2 cross section

# Sources
pulse = 0
N_photons = 5e48/(dx**3) #5e48/dx # number of photons per second
width = 1
param = 0 # photon packet or ponctual source
NS = 1


# Photo-Ionization switch
chem = 0
