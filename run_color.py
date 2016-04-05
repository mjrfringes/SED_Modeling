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
target_list = ['NGC2071.1','NGC2071.2','NGC2071.3','NGC2071.4','NGC2071.5']
#'IRAS20050.1','IRAS20050.2','IRAS20050.3','IRAS20050.4','IRAS20050.5','IRAS20050.6','IRAS20050.7',
def distance(target):
	if "IRAS20050" in target:
		return 700.*pc
	elif "NGC2071" in target:
		return 390.*pc

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

# Process data points
TwoMASS = ['j','h','ks']
uTwoMASS = ["e_"+col for col in TwoMASS]
wlTwoMASS = [1.3,1.6,2.2]
labelTwoMASS = '2MASS'
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
#d = SphericalDust()
#d.read('d03_5.5_3.0_A.hdf5')
#chi = d.optical_properties.chi
#chi = chi[::-1]
#wav = d.optical_properties.wav
#wav = wav[::-1]
#Chi = interp1d(wav,chi,kind='linear')
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

# set the wavelengths
names = TwoMASS+Spitzer+SOFIA+Herschel
wl=wlTwoMASS+wlSpitzer+wlSOFIA+wlHerschel
wl_table = Table(names = names,dtype=['f8' for col in wl])
wl_table.add_row(wl)
N = len(wl)

# load the grid
name = ['model']
folder = ['Gridfinal/']
filename = folder[0]+name[0]+'.grid.dat'
if os.path.exists(filename):
	f = open(filename,'r')
	grid = pickle.load(f)
	f.close()

else:
	#TODO: do something when the file isn't found!
	print "No grid reference file found!"
	
oldparams = ['name', 'folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
				'cav_theta','innerdustfile','outerdustfile','beta','L_sun','env_mass','env_power']
oldtypes = ['<S30','<S30','f8','f8','f8','f8','<S30','f8','f8','f8','<S30','<S30','f8','f8','f8','f8']
newgrid = Table(names=oldparams+['ext','inc']+names,dtype=oldtypes+['f8','f8']+['f8' for val in names])

# calculate value for each band for each model, extinction value and inclination
for i in range(len(grid)):
	# load model
	fname = folder[0]+grid['name'][i]+'.rtout'
	if i%10 ==0:
		print "Model: ",fname
	if os.path.exists(fname):# and grid['env_rmax'][i]==5000.0 and grid['disk_mass'][i]==0.003:
		#print "Model found!"
		mo = ModelOutput(fname)

		# load sed from model
		sed = mo.get_sed(aperture=-1, inclination='all', distance=100.*pc,units='Jy')

		for extinction in extinctions:

			# calculate optical depth
			tau_ext1 = Chi(sed.wav)/Chi(0.550)/1.086
			tau = tau_ext1*extinction
			#print "tau,",tau
	
			# calculate extinction for all inclinations
			ext = np.array([np.exp(-tau) for shape in range(sed.val.shape[0])])
			#print "ext,",ext

			# apply extinction to model
			extinct_values = np.log10(sed.val.transpose()*ext.T)
			#print "extinct_values,extinct_values.shape",extinct_values,extinct_values.shape

			# now calculate chi squared for each inclination
			for j in range(len(inclinations)):
				inc = inclinations[j]
				#print "inc=",inc
		
				# interpolate the SED from the model
				interp_func = interp1d(sed.wav,extinct_values[:,j],kind='linear')
		
				# save information into grid table
				# create a row with old parameter ples new ones
				newline = [grid[i][param] for param in oldparams]+[extinction,inc] + [interp_func(val) for val in wl]
				#print newline
		
				# add the row to the new grid
				newgrid.add_row(newline)
			
pickle.dump(newgrid,open(folder[0]+name[0]+"_eval.grid.dat",'wb'))
newgrid.write(folder[0]+name[0]+"_eval.grid.txt",format='ascii.fixed_width')



