#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::			Chemistry			::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Module importation ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
import numpy as np
import scipy.constants as sc
import config as cf
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  Constants :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c = cf.c # speed of light

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Reaction_rates ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

""" COMPILATION of reaction rates for photo-chemistry of hydrogen
taken from Maselli et al. 2003 (CRASH)"""

def gammaH0(T):
	""" DOCUMENT
	collisional ionisation rate from Maselli et al. 2003
	in cm^3/s"""
	res=5.85e-11
	res*=np.sqrt(T)
	res*=1./(1.+np.sqrt(T/1.e5))
	res*=np.exp(-157809.1/T)
	return res

def alphaAH(T):
	"""DOCUMENT
	case A recombination rate from Hui & Gnedin 1997
	in cm^3/s"""
	lambda_=2.*157807./T
	res=1.269e-13
	res*=pow(lambda_,1.503)
	res/=pow(1.+(lambda_/0.522)**0.47,1.923)
	return res

def recombination_cooling_rate_A_H(T):
	""" DOCUMENT
	case A HII recombinatino cooling rate from Hui & Gnedin 1997
	in erg cm^3 s^-1"""

	lambda_=2.*157807./T
	res=1.778e-29*pow(lambda_, 1.965)
	res/=pow(1.+pow(lambda_/0.541,0.502),2.697)
	return res

def alphaBH(T):
	""" DOCUMENT
	case B recombination rate from Hui & Gnedin 1997
	in cm^3/s
	checked against table 1 case B of Ferland et al. 1992"""

	lambda_=2.*157807./T
	res=2.753e-14
	res*=pow(lambda_,1.5)
	res/=pow(1.+(lambda_/2.74)**0.407,2.242)
	return res


def betaH(T):
	""" DOCUMENT
	HI collisional ionisation coefficient from Hui & Gnedin 1997
	in cm^3 s^-1 K^3/2 """

	lambda_=2.*157807./T
	res=21.11*pow(T,-3./2.)*np.exp(-lambda_/2.)*pow(lambda_,-1.089)
	res/=pow(1.+pow(lambda_/0.354,0.874),1.01)
	return res

def ksiH0(T):
	"""DOCUMENT
	collisional ionisation cooling from Maselli et al. 2003
	in erg cm^3 s^-1"""
	
	res=1.27e-21*np.sqrt(T)/(1.+pow(T/1.e5,0.5))
	res*=np.exp(-157809.1/T)
	return res


def etaH0(T):
	"""DOCUMENT
	recombination cooling for H0 (Maselli et al. 2003)
	case A or B or total ?
	in erg cm^3 s^-1"""

	res=8.7e-27*np.sqrt(T)*pow(T/1.e3,-0.2)/(1.+pow(T/1.e6,0.7))
	return res

def psiH0(T):
	""" DOCUMENT
	collisinal exciation cooling for H0 (Maselli et al. 2003)
	in erg cm^3 s^-1"""

	res=7.5e-19/(1.+pow(T/1.e5,0.5))
	res*=np.exp(-118348./T)
	return res

def betabremsstrahlung(T):
	"""DOCUMENT
	Bremsstrahlung cooling from Maselli et al. 2003
	in erg cm^3 s^-1
	CAREFUL: we took the densities out of the formula
	so one needs to multiply the result by rho_electrons^2 (in case of pure Hydrogen chemistry"""
	
	res=1.42e-27*np.sqrt(T)
	return res

def coolingrate(T,x):
	""" DOCUMENT
	just the sum of the cooling terms we have implemented
	just multyply by rho^2
	in erg cm^3 s^-1"""

	res=betabremsstrahlung(T)*x**2+psiH0(T)*(1.-x)**2+ksiH0(T)*(1.-x)**2+etaH0(T)*x**2
	return res
	

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Ionization functions ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def Q_root(T, N_gamma, rho, x, dt):
	""" Function who returns the value of the ionization rate at dt+1
	param: T, N_gamma, rho, x, dt
	return: X"""
	
	sigma = cf.sigma
	
	# Coeff calculation
	m = (alphaBH(T)+betaH(T))*rho*rho*dt
	n = rho - (alphaAH(T)+betaH(T))*rho/(sigma*c*1e2) -alphaBH(T)*rho*rho*dt - 2*betaH(T)*rho*rho*dt
	p = -rho*(1+x) - N_gamma -1/(sigma*c*1e2*dt) + betaH(T)*rho/(sigma*c*1e2) + betaH(T)*rho*rho*dt
	q = N_gamma + rho*x + x/(sigma*c*1e2*dt)
	
	coeff = [m,n,p,q]
	res = np.roots(coeff)
	for i in res:
		if 1>i>0:
			return i
	print("error : no root found")
	return 0

def updt_NGamma(T, N_Gamma, rho, x, dt):
	""" Function who updates N_Gamma with the chemistry module
	param: T, N_gamma, rho, x, dt
	return: N_Gamma at dt+1"""
	X = Q_root(T, N_Gamma, rho, x, dt)
	res = N_Gamma + betaH(T)*rho*rho*(1-X)*X*dt-alphaBH(T)*rho*rho*X*X*dt-rho*(X-x)
	return res*1e6, X


#::::::::::::::::::::::::::Test Chemistry ::::::::::::::::::::::::::::
#  T, N_Gamma, rho, x, dt
"""data = [[2.e3, 1.,     1.,    1.e-4, 1.e13],
				   [2.e3, 1.,     1.,    1.e-1, 1.e13],
				   [2.e3, 1.e-10, 1.,    1.e-4, 1.e13],
				   [2.e4, 1.e-10, 1.,    1.e-4, 1.e13],
				   [2.e8, 1.e-10, 1.,    1.e-4, 1.e13],
				   [2.e4, 1.,     1.,    1.e-4, 1.e13],
				   [2.e4, 1.e-2,  1.e-3, 1.e-4, 1.e13]]


print("x values :\n")

for i in range(7):
	print(Q_root(data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]))
	
print("N_Gamma values :\n")

for i in range(7):
	print(updt_NGamma(data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]))"""
	
# valeurs similaires juste diffrence pour T élevé

































