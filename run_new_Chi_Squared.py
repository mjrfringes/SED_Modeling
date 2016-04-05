'''
This script uses a simple Chi Squared estimate to estimate which model from the grid best matches the data from an SED
'''
import photometry as p
import pickle
import numpy as np
from astropy.table import Table
import os
from scipy.interpolate import interp1d
from hyperion.model import ModelOutput
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc
import matplotlib.pyplot as plt
import sys

# set the target
target_list = ['IRAS20050.1','IRAS20050.2','IRAS20050.3','IRAS20050.4','IRAS20050.5','IRAS20050.6','IRAS20050.7','NGC2071.1','NGC2071.2','NGC2071.3','NGC2071.4','NGC2071.5']
#
def distance(target):
	if "IRAS20050" in target:
		return 700.*pc
	elif "NGC2071" in target:
		return 490.*pc
		
def distpc(target):
	if "IRAS20050" in target:
		return 700.
	elif "NGC2071" in target:
		return 490.

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

# Process data points
Spitzer = ['i1','i2','i3','i4']
uSpitzer = ["e_"+col for col in Spitzer]
wlSpitzer = [3.6,4.5,5.8,8.]
labelSpitzer = 'Spitzer'
SOFIA = ['F11','F19','F31','F37']
uSOFIA = ["e_"+col for col in SOFIA]
wlSOFIA = [11.1,19.7,31.5,37.1]
labelSOFIA = 'SOFIA'
#HERSCHEL = plots(['H70','H160','H250','H350','H500'],[70,160,250,350,500],colors[7],markers[7],'HERSCHEL')
Herschel = ['H70','H160','H250','H350']
wlHerschel = [70,160,250,350]
uHerschel = ["e_"+col for col in Herschel]
labelSpitzer = 'Herschel'
sources = sourcetable.group_by('SOFIA_name')

# set up extinction
extinctions = range(30)
d = SphericalDust()
d.read('OH5.hdf5')
chi = d.optical_properties.chi#/100. # divide by 100 for the gas-to-dust ratio
chi = chi[::-1]# divide by 100 for the gas-to-dust ratio
wav = d.optical_properties.wav
wav = wav[::-1]
Chi = interp1d(wav,chi,kind='linear')

#inclinations = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]
angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
inclinations=angles[::-1]

# load the grid
name = ['model']
folder = ['Gridfinal/']
filename = folder[0]+name[0]+"_eval.grid.dat"
if os.path.exists(filename):
	f = open(filename,'r')
	grid = pickle.load(f)
	f.close()

else:
	#TODO: do something when the file isn't found!
	print "No grid reference file found!"


for target in target_list:
	print "Working on target ",target
	for key,source in zip(sources.groups.keys,sources.groups):
		if target == source['SOFIA_name'][0]:	
			# adapts to that one target and removes the J-band measurement
			if target == 'IRAS20050.1':
				TwoMASS = ['h','ks']
				uTwoMASS = ["e_"+col for col in TwoMASS]
				wlTwoMASS = [1.6,2.2]
				labelTwoMASS = '2MASS'
			else:
				TwoMASS = ['j','h','ks']
				uTwoMASS = ["e_"+col for col in TwoMASS]
				wlTwoMASS = [1.3,1.6,2.2]
				labelTwoMASS = '2MASS'

			# if the target is NGC2071, include the Herschel fluxes
			if "NGC2071" in target:
				names = TwoMASS+Spitzer+SOFIA#+Herschel
				unames = uTwoMASS+uSpitzer+uSOFIA#+uHerschel
				datapoints = source[names]
				dataerrors = source[unames]
				dataflags = source[["flag_"+col for col in names]]
				wl=wlTwoMASS+wlSpitzer+wlSOFIA#+wlHerschel
			else:
				names = TwoMASS+Spitzer+SOFIA
				unames = uTwoMASS+uSpitzer+uSOFIA
				datapoints = source[names]
				dataerrors = source[unames]
				dataflags = source[["flag_"+col for col in names]]
				wl=wlTwoMASS+wlSpitzer+wlSOFIA
		
			# calculate log10 of quantities required for chi squared
			logFnu_tot = np.log10(p.nptable(datapoints))-0.5*(1./np.log(10.))*p.nptable(dataerrors)**2/p.nptable(datapoints)**2
			varlogFnu_tot = (1./np.log(10)/p.nptable(datapoints))**2*p.nptable(dataerrors)**2
			#print "p.nptable(datapoints),p.nptable(dataerrors):",p.nptable(datapoints),p.nptable(dataerrors)


	# set the wavelengths
	wl_table = Table(names = names,dtype=['f8' for col in wl])
	wl_table.add_row(wl)
	masks = Table(names = names,dtype=['bool' for col in wl])
	# mask is 1 when flag indicates upper limit, 0 otherwise
	masks.add_row([dataflags["flag_"+col][0] =="U" for col in names])
	#print [dataflags["flag_"+col][0] =="U" for col in TwoMASS+Spitzer+SOFIA]
	#print masks
	N = len(wl)
	#fig,ax = plt.subplots()
	#ax.errorbar(wl,p.nptable(datapoints).T,p.nptable(dataerrors).T)
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	#plt.show()


		
	oldparams = ['name', 'folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
					'cav_theta','innerdustfile','outerdustfile','beta','L_sun','env_mass','env_power']
	oldtypes = ['<S30','<S30','f8','f8','f8','f8','<S30','f8','f8','f8','<S30','<S30','f8','f8','f8','f8']
	newgrid = Table(names=oldparams+['ext','inc','chi2','chi2_old','n'],dtype=oldtypes+['f8','f8','f8','f8','i8'])

	# calculate chi squared metric for each run of the grid
	#for i in range(len(grid)):
	#for i in range(1):
	# load model
	#fname = folder[0]+grid['name'][i]+'.rtout'
#	if i%10 ==0:
#		print "i: ",i
	# this calculation of Chi Squared takes into accound the upper flux limits flagged "U" in the data table
	chi2=0.0
	grid["chi2_"+target]= 0.0
	grid["n_"+target] = float(N)
	for wavel in names:
		logFnu = np.log10(datapoints[wavel][0])-0.5*(1./np.log(10.))*dataerrors["e_"+wavel][0]**2/datapoints[wavel][0]**2
		varlogFnu = (1./np.log(10)/datapoints[wavel][0])**2*dataerrors["e_"+wavel][0]**2
		intwavel = grid[wavel]-2*np.log10(distpc(target)/100.)
		# if data not upper limit, include normally in the chi2 calculation
		# when to add chi square? when not upper limit or, if upper limit, when logFnu-intwavel is negative
		grid["tmp"] = np.sign(intwavel - logFnu  )
		# when the sign is positive, then we want to add; otherwise, we don't.
		grid["tmp"][grid["tmp"]<0.0] = 0.0
		
		if ~masks[wavel][0]:
			grid["chi2_"+target] += (logFnu - intwavel)**2/varlogFnu
		else:
			# adds only when the grid"tmp" column is positive
			grid["chi2_"+target] += (logFnu - intwavel)**2/varlogFnu * grid["tmp"]
		#grid["chi2_"+target].more()
			# update the value of n
			grid["n_"+target] -= (1.0-grid["tmp"])
#			#print "chi2 not upper limit (wavel, contrib. to chi2):",wavel,(logFnu - intwavel)**2/varlogFnu
#		# if data is upper limit, then only contribute to chi2 if fit is higher than upper limit
#		else:
#			if intwavel > logFnu:
#				grid["chi2_"+target] += (logFnu - intwavel)**2/varlogFnu
#				#print "chi2 upper limit but model above measured (wavel, contrib. to chi2):",wavel,(logFnu - intwavel)**2/varlogFnu
#			else:
#				n-=1
		#grid.more()
	#print "n=",n
	grid["chi2_"+target] /= grid["n_"+target]
	#grid.sort(["chi2_"+target])
	#grid.more()
	# save information into grid table
	# create a row with old parameter ples new ones
	#newline = [grid[i][param] for param in oldparams]+[extinction,inc,chi2,chi2_old,n]
	#print newline

	# add the row to the new grid
	#newgrid.add_row(newline)
			
			#plt.loglog(sed.wav,sed.val.transpose()*ext.T, color='black')
			#plt.loglog(wl,p.nptable(datapoints).T,color='red')
	
	#plt.show()
grid.sort('chi2_NGC2071.1')
pickle.dump(grid,open(folder[0]+name[0]+"_all_chi2.grid.dat",'wb'))
grid.write(folder[0]+name[0]+"_all_chi2.grid.txt",format='ascii.fixed_width')
#newgrid.more()


