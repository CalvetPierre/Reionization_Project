#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 			Main.py		::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: Module importation ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
import chemistry as ch
import numerical_integration as num 
import matplotlib.pyplot as plt
import config as cf
import numpy as np

#::::::::::::::::::::::::::::::::::::::::      User Interface     ::::::::::::::::::::::::::::::::::::::::::::::

list_anss = ["1", "2", "3"]
list_y = ["Y", "y", "N", "n"]

def print_menu():
	print("::::::::::::::::::: Paremeters choices for the simulation ::::::::::::::::::::::\n")
	print("Chemistry:")
	ans = input("Do you want to add the chemistry module ? [y/n]")
	if ans not in list_y:
			print("error")
			return(0)
	if ans == 'y' or ans == 'Y':
		cf.chem = 1
		print("The chemistry module implemented")
	else:
		print("No chemistry module")
	print("\n")
	print("Choice of Sources: ")
	print("1. Pulse ")
	print("2. Continuous sources ")
	print("3. Packet source ")
	anss = input("Please enter your choice: ")
	print(anss)
	if anss not in list_anss:
		print("error")
	else: 
		if anss == "1":
			cf.pulse = 1 # les autres paramètres par défaut sont ok
		if anss == "2":
			ansm = input("Do you want multiple sources ? [y/n]")
			if ansm not in list_y:
				print("error")
				return(0)
			if ansm == "y" or ansm == "Y":
				ansn = input("Please enter the number of multiple sources you want :")
				if ansn.isdigit():
					print("ok")
				cf.NS = int(ansn)
				print("You want", str(cf.NS) , " sources")
		if anss == "3" :
			cf.param = 1
			answ = input("Please enter the width of your source packet in term of cells: ")
			cf.width = int(answ)
			print("You want a width of ", str(cf.width) , " cells")
	return(0)
	
#:::::::::::::::::::::::::::::::::::::::: Test Numerical Integration :::::::::::::::::::::::::::::::::::::::::::

print_menu()

print("Theoretical Stromgren length:", cf.N_photons*cf.dx/(ch.alphaBH(cf.T)*cf.rho*cf.rho*1e6))
print("Length:", cf.L)
print("Total distance reachable by photon:",cf.tf*cf.c)
print("Source emission:", cf.N_photons)
print("Total time:", cf.tf)
print("Temperature:", cf.T)
print("rho:",cf.rho)
print("Chimie:",cf.chem)



for i in range(cf.N_t):
	num.GLF()
	num.source(cf.param, cf.NS, cf.width)
	num.updt_N(cf.chem)
	num.updt_F(cf.chem)
	num.compute_khi()
	num.compute_P()
	"""num.compute_T()"""
	if (i%(cf.N_t/10) == 0):
		print(i,"/",cf.N_t)
		plt.figure(1)
		plt.title("Photon density for a ponctual source (with Hydrogen photo-ionization")
		plt.plot(np.arange(0, cf.N_c, 1)*cf.dx, cf.N)
		# plt.axvline((cf.N_c/2+6)*cf.dx + (i*cf.dt*cf.c))
		# plt.axvline(cf.N_c/2*cf.dx + (cf.N_photons*cf.dx/(ch.alphaBH(cf.T)*cf.rho*cf.rho*1e6)/2))
		# plt.axvline(cf.N_c/2*cf.dx + cf.tf*cf.c)
		plt.xlabel(" Distance x [m]")
		plt.ylabel("Photon Density [/m³]")
		
		
		if cf.chem == 1:
			plt.figure(3)
			plt.title("Ionization rate")
			plt.plot(np.arange(0, cf.N_c, 1)*cf.dx, cf.x)
			# plt.axvline(cf.N_c/2*cf.dx + cf.N_photons*cf.dx/(ch.alphaBH(cf.T)*cf.rho*cf.rho*1e6))
			# plt.axvline(cf.N_c/2*cf.dx + cf.tf*cf.c)

plt.show()

"""
plt.figure()
plt.plot(cf.N)
plt.show()

plt.figure()
plt.plot(cf.P, label = "Pressure", marker="o")
plt.plot(cf.P_GLF, label = "Pressure GLF", marker="o")
plt.legend()
plt.show()

plt.figure()
plt.plot(cf.F, label="Flux", marker="o")
plt.plot(cf.F_GLF, label="Flux GLF",marker="o")
plt.legend()
plt.show()

plt.figure()
plt.plot(cf.Khi)
plt.show()
"""

print("\n End!")
