#! /usr/bin/env python3

import numpy as np
import math
import os

from spectrum_2v2b import me,Qbb,keV,MeV,M2nuEff
from spectrum_2v2b import dG0dEe1dEe2,dG0dEe1,dG1dEe1,dWdT

# Set neutrino mass values
mnu1 = 400*keV
mnu2 = 400*keV

# Do you only neeed the summed spectrum?
dosumonly = False
# Do you want to compute 2D grid? it will take a while
do2D = False

## Choose your mother nucleus
mother = '76Ge'
#mother = '136Xe'
#mother = '100Mo'

print('\nComputing spectra for '+ mother + ' with mnu1='
	+str(int(mnu1/keV))+'keV and mnu2=' +str(int(mnu2/keV)) +'keV')
print('\nComputing 1D spectra. It might take a couple of minutes...')
if (dosumonly):
	print('Computing only summed spectrum. For more results, set dosumonly = False')


# where to export the results
pathresult   = mother+'_2vbb_mnu_'+str(int(mnu1))+'_'+str(int(mnu2))+'_keV/'

if not os.path.exists(pathresult): 
	os.mkdir(pathresult)

if (not dosumonly):
	fses  = open(pathresult+mother+'_ses.txt','w')
	fcor  = open(pathresult+mother+'_cor.txt','w')
fsums = open(pathresult+mother+'_sums.txt','w')

# running the loop
# T = E-m kinetic energy
Tmin 	= 1
Tmax 	= Qbb
Tpoints = int(Qbb)
Tlist	= np.linspace(Tmin,Tmax,Tpoints)
i=0 #counter
for T1 in Tlist:
	i+=1

	# single electron spectrum
	dG0   = dG0dEe1(me+T1,mnu1,mnu2)
	dWdT1 = M2nuEff*M2nuEff*dG0

	if (not dosumonly):
	# angular correlation
		if (dG0==0): # Set to zero when no PS available, as in 1209.5722
			alphaE = 0
		else:
			dG1    = dG1dEe1(me+T1,mnu1,mnu2)
			alphaE = dG1/dG0

	# summed energy spectrum
	dWdT12 = dWdT(T1,mnu1,mnu2) # here T1 means T=T1+T2

	# export output
	if (not dosumonly):
		fses.write(('%4d %.4f %.4e \n') % (i,T1/MeV,dWdT1))
		fcor.write(('%4d %.4f %.4e \n') % (i,T1/MeV,alphaE))
	fsums.write(('%4d %.4f %.4e \n') % (i,T1/MeV,dWdT12))

print('1D results exported to '+ pathresult)
if (not dosumonly):
	fses.close()
	fcor.close()
fsums.close()


#2D grid
if (do2D and not dosumonly):
	print('\n Computing 2D spectra. It might take a bit longer...')

	# T = E-m kinetic energy
	Tmin 	= 1
	Tmax 	= Qbb-1
	Tpoints = int(Qbb)
	Tlist	= np.linspace(Tmin,Tmax,Tpoints)

	fdata = open(pathresult+mother+'_2ds.txt','w')
	i1=0
	for T1 in Tlist:
		i1+=1
		i2=0
		for T2 in Tlist:
			i2+=1
			if(T1+T2>=Qbb-mnu1-mnu2):
				fdata.write(('%d %4d %.4f %.4f %.4e \n') % (i1,i2,T1/MeV,T2/MeV,0))
				break
			else:
				dW = dG0dEe1dEe2(me+T1,me+T2,mnu1,mnu2)
				fdata.write(('%d %4d %.4f %.4f %.4e \n') % (i1,i2,T1/MeV,T2/MeV,dW))
	    	
    
	fdata.close()
	print('2D results exported to '+ pathresult)

else:
	print('\n Skipping 2D spectra. Switch do2d=True in get_spectrum_2v2b.py to compute it.')




print('done\n')
