#:::::::::::::::::::::::::::::::::::::::::::::::::: 			Numerical Integration		::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Module importation ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
import chemistry as ch
import scipy.constants as sc
import config as cf
import numpy as np
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  Constants :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c = cf.c # speed of light

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Integration module ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

""" Steps order : 
GLF P and F because m = n (no specific order) 
Update N and F for n+1 (no specific order) 
Compute khi for n+1
Compute P for n+1"""

# GLF Function

def GLF():
	""" function computes the GLF for P and F
	params: none
	return : none
	GLF Construction:
	Classic:  #i-1   #i   #i+1
	GLF    :#i-1   #i   # i+1"""
	for i in range(cf.N_c):
		cf.F_GLF[i] = 0.5*((cf.F[i-1]+cf.F[i])-c*(cf.N[i]-cf.N[i-1]))
		cf.P_GLF[i] = 0.5*((cf.P[i-1]+cf.P[i])-c*(cf.F[i]-cf.F[i-1]))
		
# Update Functions

def updt_N(chem = 0):
	""" function updates N at t = n+1
	param : chem if chem = 0 : no chemistry if chem = 1 : chemistry module
	return: none"""
	for i in range(cf.N_c-1):
		cf.N[i] += (cf.dt/cf.dx)*(cf.F_GLF[i]-cf.F_GLF[i+1])
	cf.N[cf.N_c-1] += (cf.dt/cf.dx)*(cf.F_GLF[cf.N_c-1]-cf.F_GLF[0])
	
	if chem == 1: # chemistry module taken into account
		for i in range(cf.N_c):
			cf.N[i], cf.x[i] = ch.updt_NGamma(cf.T, cf.N[i]*1e-6, cf.rho, cf.x[i], cf.dt)

def updt_F(chem = 0):
	""" function updates F at t = n+1
	param : none
	return: none"""
	for i in range(cf.N_c-1):
		cf.F[i] += (cf.dt/cf.dx)*(cf.P_GLF[i]-cf.P_GLF[i+1])
	cf.F[cf.N_c-1] += (cf.dt/cf.dx)*(cf.P_GLF[cf.N_c-1]-cf.P_GLF[0])
	
	if chem == 1: # chemistry module taken into account
		for i in range(cf.N_c):
			cf.F[i] /= (1+ cf.c*cf.sigma*1e4*cf.rho*1e6*cf.dt*(1-cf.x[i]))

# Compute Functions

def compute_khi():
	""" function computes khi
	param: none
	return: none"""
	
	f = cf.F/(c*cf.N)
	for i in range(len(f)):
		if (f[i] > np.sqrt(4/3)):
			print(f[i])
			f[i] = np.sqrt(4/3)
		if (f[i] < -np.sqrt(4/3)):
			print(f[i])
			f[i] = -np.sqrt(4/3)
	cf.Khi = (3+4*f*f)/(5+2*np.sqrt(4-3*f*f))

def compute_P():
	""" function computes P
	param: none
	return: none"""
	cf.P = cf.Khi*cf.N*c*c

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Sources ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def updt_NSources(NS):
	""" function that add the photon source term
	param: NS the number of sources
	return: none"""
	for i in range(NS):
		cf.N[int(cf.N_c/(1+NS))*(i+1)] += cf.N_photons*(cf.dt)
		
	if cf.pulse == 1:
		cf.N_photons=0

def updt_PacketSource(width):
	""" function that add the photon source term
	param: width the width of the source in term of cells number
	return: none"""
	cf.N[int(cf.N_c/2-width/2):int(cf.N_c/2+width/2)] += cf.N_photons*(cf.dt)

def source(param = 0, NS = 1, width = 1):
	if param == 0:
		updt_NSources(NS)
	else :
		updt_PacketSource(width)
		
		
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Temperature coupling ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
"""		
def compute_T():
	
	cf.T = cf.T/2 + (cf.H-cf.L)*cf.dt/(3*cf.rho*(1+cf.x)*cf.kb)
"""	
	
	
	
	
		
