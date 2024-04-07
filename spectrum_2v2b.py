import numpy as np 
import math
from scipy.integrate import quad
from scipy.special import gamma 
import sys

# Define units
keV   = 1. 					# set to 1 to work in keV
MeV   = 1e3*keV
GeV   = 1e6*keV

year  = 3.16e7/6.58e-22/MeV # 1 year in MeV-1
fm = 1./197/MeV; 			# fm in MeV-1
pi    = np.pi

# Define parameters
me    = 511.*keV 		# electron mass. PDG'20
gA    = 1.27 			# axial coupling. PDG'20
Gb    = 1.09e-5/GeV/GeV # Gfermi * cos\theta_Cabibbo. PDG'20
alpha = 1./137.			# alpha_EM


# mother to be set at get_spectrum_2v2b.py main file

if (mother == '76Ge'):
	Qbb     = 2039*keV		# End-point 
	Z       = 32. 			# Atomic number 
	A 		= 76 			# Mass number
	Atilde  = 9.411*MeV 	# Atilde parameter [1303.4124]
	M2nuEff = 0.118 		# eff NME from Table VI of 1209.5722
	R 		= A**(1./3)*1.2*fm
	#R 	= 2.66e-5/keV 	     # Using R=R(A) of Nuclear Physics 5, 173 (1958)

elif (mother == '136Xe'):
	Qbb     = 2458*keV		# End-point 
	Z		= 54. 			# Atomic number 
	A 		= 136 			# Mass number
	Atilde  = 13.06*MeV 	# Atilde parameter [1303.4124]
	M2nuEff = 0.0182 		# eff NME from Table VI of 1209.5722
	R 		= A**(1./3)*1.2*fm
	#R 	= 3.14e-5/keV 		# Using R=R(A) of Nuclear Physics 5, 173 (1958)

elif (mother == '100Mo'):
	Qbb     = 3034*keV		# End-point 
	Z		= 42. 			# Atomic number 
	A 		= 100 			# Mass number
	Atilde  = 11.20*MeV 	# Atilde parameter [1303.4124]
	M2nuEff = 0.206 		# eff NME from Table VI of 1209.5722
	R 		= A**(1./3)*1.2*fm
	#R 	= 2.87e-5/keV 		# Using R=R(A) of Nuclear Physics 5, 173 (1958)

else:
	print('\nERROR 404')
	print('Mother nucleus ' + mother + ' not found. Please check spectrum_2v2b.py.')
	print('Program will Exit.\n')
	sys.exit()

# R=R(A) as in Nuclear Physics 5, 173 (1958)
def nuclR(A): 
	fm = 1./197/MeV;
	return (1.121*A**(1./3) + 2.426*A**(-1./3) - 6.6144/A)*fm;

# Phase space including neutrino masses
def w2nu(Ee1,Ee2,Enu1,mnu1,mnu2):
	# apply energy conservation for Enu2
	Enu2 = Qbb+2*me-Ee1-Ee2-Enu1

	# define momenta
	pe1  = np.sqrt(Ee1*Ee1-me*me)
	pe2  = np.sqrt(Ee2*Ee2-me*me)
	pnu1 = np.sqrt(Enu1*Enu1-mnu1*mnu1)
	pnu2 = np.sqrt(Enu2*Enu2-mnu2*mnu2)

	return gA**4*Gb**4/64/pi**7*pe1*Ee1*pe2*Ee2*pnu1*Enu1*pnu2*Enu2

# Fermi function
def F0(E): 
	Zf    = Z+2 # Atomic number of final daugther nucleus 
	p     = np.sqrt(E*E-me*me)
	gam   = np.sqrt(1.-alpha*alpha*Zf*Zf) # gamma parameter, not the same as gamma()=Euler Gamma function
	y	  = alpha*Zf*E/p

	#return 2*pi*y/(1-np.exp(-2*pi*alpha*Zf))   # Primakoff-Rosen approx		
	#return 2*pi*y/(1-np.exp(-2*pi*y))			# Non-relativistic approx
	return 4*(2*p*R)**(2*gam-2)*np.exp(pi*y)*(np.abs(gamma(gam+y*1j))/gamma(2*gam+1))**2 # Full relativistic expression

def f110(Ee1,Ee2):
	return F0(Ee1)*F0(Ee2)

def f111(Ee1,Ee2):
	# define momenta
	pe1  = np.sqrt(Ee1*Ee1-me*me)
	pe2  = np.sqrt(Ee2*Ee2-me*me)

	return -pe1*pe2/Ee1/Ee2*F0(Ee1)*F0(Ee2)


def KN(Ee1,Ee2,Enu1):
	# apply energy conservation for Enu2
	Enu2  = Qbb+2*me-Ee1-Ee2-Enu1
	# estimate <EN>-EI = Atilde - (Qbb - 2me)/2
	Ediff = Atilde - Qbb/2 - me
	return 1./(Ee1+Enu1+Ediff) + 1./(Ee2+Enu2+Ediff)

def LN(Ee1,Ee2,Enu1):
	# apply energy conservation for Enu2
	Enu2  = Qbb+2*me-Ee1-Ee2-Enu1
	# estimate <EN>-EI = Atilde - (Qbb - 2me)/2
	Ediff = Atilde - Qbb/2 - me
	return 1./(Ee1+Enu2+Ediff) + 1./(Ee2+Enu1+Ediff)


# Full differential rate for G0
def dG0(Ee1,Ee2,Enu1,mnu1,mnu2):
	KNval = KN(Ee1,Ee2,Enu1)
	LNval = LN(Ee1,Ee2,Enu1)
	KL0   = KNval*KNval + LNval*LNval + KNval*LNval
	fact0 = 2*Atilde**2/(gA**4*me**2*3*np.log(2))*year
	return fact0*f110(Ee1,Ee2)*KL0*w2nu(Ee1,Ee2,Enu1,mnu1,mnu2)

# Full differential rate for G1
def dG1(Ee1,Ee2,Enu1,mnu1,mnu2):
	KNval = KN(Ee1,Ee2,Enu1)
	LNval = LN(Ee1,Ee2,Enu1)
	KL1   = 2*KNval*KNval + 2*LNval*LNval + 5*KNval*LNval
	fact1 = 2*Atilde**2/(gA**4*me**2*9*np.log(2))*year
	return fact1*f111(Ee1,Ee2)*KL1*w2nu(Ee1,Ee2,Enu1,mnu1,mnu2)

# Differential of G0 rate wrt Ee1 and Ee2
def dG0dEe1dEe2(Ee1,Ee2,mnu1,mnu2):
	Enu1min = mnu1
	Enu1max = Qbb+2*me-Ee1-Ee2-mnu2
	return quad(lambda Enu1: dG0(Ee1,Ee2,Enu1,mnu1,mnu2),Enu1min,Enu1max)[0]

# Differential of G1 rate wrt Ee1 and Ee2
def dG1dEe1dEe2(Ee1,Ee2,mnu1,mnu2):
	Enu1min = mnu1
	Enu1max = Qbb+2*me-Ee1-Ee2-mnu2
	return quad(lambda Enu1: dG1(Ee1,Ee2,Enu1,mnu1,mnu2),Enu1min,Enu1max)[0]

# Differential of G0 rate wrt Ee1 
def dG0dEe1(Ee1,mnu1,mnu2):
	if (Ee1>=Qbb+me-mnu1-mnu2):
		return 0
	else:
		Ee2min = me
		Ee2max = Qbb+2*me-Ee1-mnu1-mnu2
		return quad(lambda Ee2: dG0dEe1dEe2(Ee1,Ee2,mnu1,mnu2),Ee2min,Ee2max)[0]

# Differential of G1 rate wrt Ee1 
def dG1dEe1(Ee1,mnu1,mnu2):
	if (Ee1>=Qbb+me-mnu1-mnu2):
		return 0
	else:
		Ee2min = me
		Ee2max = Qbb+2*me-Ee1-mnu1-mnu2
		return quad(lambda Ee2: dG1dEe1dEe2(Ee1,Ee2,mnu1,mnu2),Ee2min,Ee2max)[0]

# Differential of W2nu wrt Ee1+Ee2-2me
def dWdT(T,mnu1,mnu2):
	# T = E1+E2-2*me the sum of kinetic energies
	if (T>=Qbb-mnu1-mnu2):
		return 0
	else:
		T1min = 0
		T1max = T
		return quad(lambda T1: M2nuEff*M2nuEff*dG0dEe1dEe2(me+T1,me+T-T1,mnu1,mnu2),T1min,T1max)[0]


# Differential of G0 rate wrt Ee1 (now w/o if, so its integral works better)
def dG0dEe1int(Ee1,mnu1,mnu2):
	Ee2min = me
	Ee2max = Qbb+2*me-Ee1-mnu1-mnu2
	return quad(lambda Ee2: dG0dEe1dEe2(Ee1,Ee2,mnu1,mnu2),Ee2min,Ee2max)[0]

# Phase space factor (integrated)
def G0(mnu1,mnu2):
	Ee1min = me
	Ee1max = Qbb+me-mnu1-mnu2
	return quad(lambda Ee1: dG0dEe1int(Ee1,mnu1,mnu2),Ee1min,Ee1max)[0]
