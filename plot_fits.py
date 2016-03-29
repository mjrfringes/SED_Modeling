'''
this little script plots a dataset and one or more fits with it
First argument is the name of the target data points to print
arguments 2-N are the various fits to plot
Usage (plots the IRAS20050.4 data points and the number 168 : 
python plot_fits.py IRAS20050.4 IRAS20050_168
'''


import matplotlib.pyplot as plt
import photometry as phot
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


target = sys.argv[1]
plotlist = sys.argv[2:]

dist = 700*pc

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

markers = ['v','p','D','^','h','o','*','>','<']
TwoMASS = phot.plots(['j','h','ks'],[1.3,1.6,2.2],sns.xkcd_rgb['pale red'],markers[0],'2MASS')
Spitzer = phot.plots(['i1','i2','i3','i4','m1','m2'],[3.6,4.5,5.8,8.,24,70],sns.xkcd_rgb['medium green'],markers[1],'Spitzer')
SOFIA = phot.plots(['F11','F19','F31','F37'],[11.1,19.7,31.5,37.1],sns.xkcd_rgb['denim blue'],markers[3],'SOFIA')
fluxlist = [TwoMASS,Spitzer,SOFIA]
fluxnames = [p.bands for p in fluxlist]

fig,ax = plt.subplots()

sources = sourcetable.group_by('SOFIA_name')
for key,source in zip(sources.groups.keys,sources.groups):
	if target == source['SOFIA_name'][0]:	
		#datapoints = source[TwoMASS+Spitzer+SOFIA]
		#dataerrors = source[uTwoMASS+uSpitzer+uSOFIA]
		for fl in fluxlist:
			phot.plotData(ax,source,fl,0.8,msize=12)


#ax.errorbar(wl,p.nptable(datapoints).T,p.nptable(dataerrors).T,linetype=None)
ax.set_xscale('log')
ax.set_yscale('log')

# set up extinction
chi = np.loadtxt('kmh94_3.1_full.chi')
wav = np.loadtxt('kmh94_3.1_full.wav')
Chi = interp1d(wav,chi,kind='linear')

# now load up the grid
name = ['IRAS20050']
folder = ['/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/']
filename = folder[0]+name[0]+"_"+target+".grid.dat"
if os.path.exists(filename):
	grid = pickle.load(open(filename,'r'))
	# for each fit that is desired, search for the right line in the grid
	fitnames = grid.group_by('name')
	for fitkey,fitname in zip(fitnames.groups.keys,fitnames.groups):
		name = fitkey['name']
		#print name
		fname = folder[0]+name+'.rtout'
		mo = ModelOutput(fname)
		sed = mo.get_sed(aperture=-1, inclination='all', distance=dist,units='Jy')
		
		if name in plotlist:
			# sort table according to chi2
			fitname.sort('chi2')
			
			# the first line is then the best fit. let's select the extinction
			extinction = fitname['ext'][0]
			print "extinction = ",extinction
			
			# inclination
			angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
			incs=angles[::-1]
			#incs = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]

			inc = int(np.argwhere(incs==fitname['inc'][0]))
			print inc
			
			# calculate optical depth
			tau_ext1 = Chi(sed.wav)/Chi(0.550)/1.086
			print "sed.wav = ",sed.wav
			print "tau_ext1 = ",tau_ext1
			tau = tau_ext1*extinction
			print "tau,",tau
			
			# get the sed
			sed = mo.get_sed(aperture=-1, inclination=inc, distance=dist)
			
			# calculate extinction for all inclinations
			ext = np.array(np.exp(-tau))
			ax.plot(sed.wav,sed.val*ext.T, color='red',lw=3)
			#ext = np.array([np.exp(-tau) for shape in range(sed.val.shape[0])])
			#ax.plot(sed.wav,sed.val.T*ext.T, color='red',lw=3)
			
			
		
else:
	print "No grid reference file found!"





plt.show()

