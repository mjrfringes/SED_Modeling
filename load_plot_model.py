''' this loads up a pre-run model and plots things about it '''


import models as m
import pickle
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
from hyperion.model import ModelOutput

#print sys.argv
#model_folder = sys.argv[1]
#grid_file = sys.argv[2]
#grid = pickle.load(open(model_folder+grid_file,'r'))
#model_name = sys.argv[3:]

model_folder = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/"
grid_file = "model.grid.dat"
f = open(model_folder+grid_file,'r')
grid = pickle.load(f)
f.close()

model_name = ["model_402"]


amb_dens = 1e-20
inc=0
from hyperion.model import ModelOutput
angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
angles=angles[::-1]
print grid
if len(model_name)<2:
	for i in range(len(grid)):
		if grid['name'][i] == model_name[0]:
			print "test"
			model = m.YSOModelSim(name=grid['name'][i],folder=grid['folder'][i],L_sun=grid['L_sun'][i],M_sun=grid['M_sun'][i],env_mass=grid['env_mass'][i],\
				env_type='power',disk_rmax=grid['disk_rmax'][i],env_rmin=grid['env_rmin'][i],env_power=grid['env_power'][i],\
				disk_mass=grid['disk_mass'][i],cav=True,cav_theta=grid['cav_theta'][i],\
				amb_rmin=1,amb_rmax=grid['env_rmax'][i],env_rmax=grid['env_rmax'][i],T=grid['T'][i],innerdustfile=grid['innerdustfile'][i],\
				angles=angles,angles2=[0. for a in angles],env=True,disk=grid['disk'][i],\
				outerdustfile=grid['outerdustfile'][i])
			model.initModel()
			print "test"
			#model.runModel()
			model.plotSim(extinction=0,show=True,inc=inc,dist_pc=100)

#else:
#fig = plt.figure(figsize=(15,7.5))
#ax = fig.add_subplot(111)
#gs = gridspec.GridSpec(1,len(model_name))
#gs.update(wspace=0.0,hspace=0.0)
#colors = sns.color_palette('husl',9)
#for i in range(len(model_name)):
#		model = model_folder+model_name[i]
#		mo = ModelOutput(model)
#		sed = mo.get_sed(aperture=-1, inclination='all', distance=100,units='Jy')
#		print sed
#		ax.plot(sed.wav,sed.val.transpose(),color=colors[i],alpha=0.8)

#ax.set_xscale('log')
#ax.set_yscale('log')

plt.show()

