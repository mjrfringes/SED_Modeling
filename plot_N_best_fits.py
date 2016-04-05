'''
this little script plots a dataset and one or more fits with it
First argument is the name of the target data points to print
arguments 2-N are the various fits to plot
Usage (plots the IRAS20050.4 data points and the 5 best fits available) : 
python plot_fits.py 10 IRAS20050.4 IRAS20050.5
'''


import matplotlib.pyplot as plt
import photometry as phot
from matplotlib.ticker import ScalarFormatter
from hyperion.dust import SphericalDust
import pickle
import numpy as np
import sys
from hyperion.model import ModelOutput
from hyperion.util.constants import pc
import os
import seaborn as sns
import matplotlib as mpl
from scipy.interpolate import interp1d

# load plot params
sns.set(context='paper',style='ticks')
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.labelsize']=16
mpl.rcParams['legend.fontsize']=14
mpl.rcParams['font.size']=14


targetlist = sys.argv[2:]
Nbest = int(sys.argv[1])

def distance(target):
	if "IRAS20050" in target:
		return 700.*pc
	elif "NGC2071" in target:
		return 390.*pc

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

markers = ['v','p','D','^','h','o','*','>','<']
TwoMASS = phot.plots(['j','h','ks'],[1.3,1.6,2.2],sns.xkcd_rgb['pale red'],markers[0],'2MASS')
Spitzer = phot.plots(['i1','i2','i3','i4','m1','m2'],[3.6,4.5,5.8,8.,24,70],sns.xkcd_rgb['medium green'],markers[1],'Spitzer')
SOFIA = phot.plots(['F11','F19','F31','F37'],[11.1,19.7,31.5,37.1],sns.xkcd_rgb['denim blue'],markers[3],'SOFIA')
HERSCHEL = phot.plots(['H70','H160','H250','H350'],[70.,160.,250.,350.],sns.xkcd_rgb['purple'],markers[4],'HERSCHEL')
VANKEMPEN = phot.plots(['SMA1300'],[1300],sns.xkcd_rgb['denim blue'],markers[5],'VANKEMPEN')
fluxlist = [TwoMASS,Spitzer,HERSCHEL,SOFIA,VANKEMPEN]
fluxnames = [p.bands for p in fluxlist]


# set up extinction
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


# now load up the grid
name = ['model']
folder = ['/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/']
filename = folder[0]+name[0]+"_all_chi2.grid.dat"
f = open(filename,'r')
grid = pickle.load(f)
f.close()

# select columns to display in console
display = ['name', 'folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
					'cav_theta','beta','L_sun','env_mass','env_power']
sources = sourcetable.group_by('SOFIA_name')
for target in targetlist:
	print target
	fig,ax = plt.subplots()

	# labels
	ax.set_xlabel(r'Wavelength (microns)')
	ax.set_ylabel(r'$\lambda F_\lambda$ (ergs.cm$^{-2}$.s$^{-1}$)')

	#ax.errorbar(wl,p.nptable(datapoints).T,p.nptable(dataerrors).T,linetype=None)
	ax.set_xscale('log')
	ax.set_yscale('log')
	
	# necessary to not have log notation
	ax.xaxis.set_major_formatter(ScalarFormatter())
	#ax.set_xticks([1,10,100,1000])

	# list of angles used when running the models
	angles=np.arccos(np.linspace(0,1.,20))*180./np.pi

	colors=	sns.cubehelix_palette(Nbest,reverse=True)
	
	# sort grid according to chi2 for given target
	grid.sort("chi2_"+target)
	print grid[:Nbest][display+['ext','inc',"chi2_"+target,"n_"+target]]

	# loop over the number of fits to plot
	for i in reversed(range(Nbest)):

		# load model
		name = grid['name'][i]
		fname = folder[0]+name+'.rtout'
		print fname
		mo = ModelOutput(fname)
	
		# load extinction
		extinction = grid['ext'][i]
		#print "extinction = ",extinction

		# load inclination
		incs=angles[::-1]
		inc = int(np.argwhere(incs==grid['inc'][i]))
		#print inc
	
		# get the sed
		sed = mo.get_sed(aperture=-1, inclination=inc, distance=distance(target),uncertainties=True)
		#print name
		#print sed.unc/sed.val
		#plt.plot(sed.wav,sed.unc/sed.val)
		#plt.show()

		# calculate optical depth
		tau_ext1 = Chi(sed.wav)/Chi(0.550)/1.086
		#print "sed.wav = ",sed.wav
		#print "tau_ext1 = ",tau_ext1
		tau = tau_ext1*extinction
		#print "tau,",tau

		# calculate extinction for all inclinations
		ext = np.array(np.exp(-tau))
		#print "exp(tau):",ext
		if i==0:
			ax.plot(sed.wav,sed.val*ext.T, lw=5,alpha=0.9,color=sns.xkcd_rgb['amber'])
			#ax.plot(sed.wav,(sed.val+sed.unc)*ext.T,lw=5,alpha=0.9,color=sns.xkcd_rgb['amber'])
			#ax.plot(sed.wav,(sed.val-sed.unc)*ext.T,lw=5,alpha=0.9,color=sns.xkcd_rgb['amber'])
			sed2 = mo.get_sed(aperture=-1, inclination=inc, distance=distance(target),units='Jy')
			np.savetxt(folder[0]+target+"_SED.txt",np.array([sed2.wav,sed2.val*ext.T]).T)
		else:
			ax.plot(sed.wav,sed.val*ext.T, lw=1,alpha=0.6,color='gray')
#		ax.plot(sed.wav,sed.val*ext.T, lw=2,alpha=1,color=colors[i])
		#ext = np.array([np.exp(-tau) for shape in range(sed.val.shape[0])])
		#ax.plot(sed.wav,sed.val.T*ext.T, color='red',lw=3)
			
	for key,source in zip(sources.groups.keys,sources.groups):
		if target == source['SOFIA_name'][0]:	
			for fl in fluxlist:
				#`print fl.bands
				phot.plotData(ax,source,fl,0.8,msize=12)

	s = '%.2f' % (grid['chi2_'+target][0])
	fig.suptitle(target + r'   Best $\chi^2 =$' + s)
	fig.tight_layout()
	fig.savefig(folder[0]+target+"_"+str(Nbest)+".png",dpi=300)
	#plt.show()

