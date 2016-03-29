''' this function exports the spectra corresponding to the best fits for a list of targets'''

from hyperion.model import ModelOutput
from hyperion.dust import SphericalDust
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pickle
from hyperion.util.constants import pc


target_list = ['IRAS20050.1','IRAS20050.2','IRAS20050.3','IRAS20050.4','IRAS20050.5']
dist = 700*pc


folder = ['Grid/']
name = ['model']

angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
inclinations=angles[::-1]

d = SphericalDust()
d.read('d03_5.5_3.0_A.hdf5')
chi = d.optical_properties.chi
chi = chi[::-1]
wav = d.optical_properties.wav
wav = wav[::-1]
Chi = interp1d(wav,chi,kind='linear')

sorted_grid = pickle.load(open(folder[0]+name[0]+"_"+target+".grid.dat",'r'))
best_model_fname = folder[0]+sorted_grid['name'][0]+'.rtout'
best_model = ModelOutput(fname)
inc = int(np.argwhere(inclinations==sorted_grid['inc'][0]))
sed = best_model.get_sed(aperture=-1, inclination=inc, distance=dist,units='Jy')
N = len(sed.wav)
vec = np.zeros(N,len(target_list)+1)
vec[:,0] = sed.wav

for i in range(len(target_list)):
	target = target_list[i]
	sorted_grid = pickle.load(open(folder[0]+name[0]+"_"+target+".grid.dat",'r'))
	best_model_fname = folder[0]+sorted_grid['name'][0]+'.rtout'
	best_model = ModelOutput(fname)
	extinction = sorted_grid['ext'][0]
				
	# get inclination
	inc = int(np.argwhere(inclinations==sorted_grid['inc'][0]))

	# get sed for best fit
	sed = best_model.get_sed(aperture=-1, inclination=inc, distance=dist,units='Jy')
	
	tau_ext1 = Chi(sed.wav)/Chi(0.550)/1.086
	tau = tau_ext1*extinction
	#print "tau,",tau

	# calculate extinction for all inclinations
	ext = np.exp(-tau)#np.array([np.exp(-tau) for shape in range(sed.val.shape[0])])
	#print "ext,",ext

	# apply extinction to model
	extinct_values = sed.val*ext

	vec[:,i+1] = extinct_values

print vec

