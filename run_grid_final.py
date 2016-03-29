'''
This program runs hyperion on a grid or predefined models
'''
import itertools
from astropy.table import Table
import models as m
import numpy as np
import pickle,os
from hyperion.util.constants import rsun, lsun, au, msun, yr, c, pc, sigma

### these are the parameters to probe:
# make sure all quantities are entered as floating points (e.g. 1.0 instead of 1)
disk = ["Flared"]
name = ['model']
folder = ['Gridfinal/']
T_list=[5000.]
m_sun_list = [1.]
l_sun_list =  [120.]#np.arange(20.,130.,5.)
disk_mass_list = [0.2]
disk_rmax_list = [200.]
#disk_rmin_list = [1.]
disk_h0_list = [0.1]


env_power_list = [-2.0]#[-1.0,-1.5,-2.0]

env_rmax_list = [5000.]

env_rmin_list = [30] #0.272
beta_list = [1.1]#[1.,1.1,1.2,1.25]

env_rho0 = 2.52e-17 # nominal density at r0 from model 63
alpha=3.0+env_power_list[0]
Mconst = 4.0*np.pi*env_rho0/alpha/(env_rmin_list[0]*au)**env_power_list[0]/msun
#env_mass_list = [0.001,0.01,0.1,0.2,0.5,1.0,2.0,4.0,8.0,12.0]#[Mconst*((Rmax*au)**alpha - (env_rmin_list[0]*au)**alpha) for Rmax in env_rmax_list]#[0.05,0.01]#[0.2,0.5,1.,2.,5.,0.1,0.05]
env_mass_list = [15.]#[0.5*(env_rmax_list[0]/5000.)**alpha]
print env_mass_list

mdot_list = [1e-6]
cav_theta_list = [49]

innerdust_list = ['OH5.hdf5']
outerdust_list = ['OH5.hdf5']

# this creates the list of lists
list_of_lists = [name, folder, T_list,m_sun_list,env_rmax_list,env_rmin_list,disk, disk_mass_list, disk_rmax_list,cav_theta_list,  innerdust_list,outerdust_list,beta_list,l_sun_list,  env_mass_list, env_power_list]

# flatten the list of lists
list_of_lists_flattened = list(itertools.product(*list_of_lists))
print list_of_lists

# inclination angles:
angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
angles=angles[::-1]

# if file exists, load existing grid file
filename = folder[0]+name[0]+'.grid.dat'
if os.path.exists(filename):
	# for all the lists, we need to check if the list already exists in the grid file
	for mylist in list_of_lists_flattened:	
		# by default, we will add a line
		add = 'yes'
		
		# open file
		f = open(filename,'r')
		grid = pickle.load(f)
		f.close()
		print grid['name'][-1]
		
		# load up last line of grid
		Ngrid = int(grid['name'][-1].split('_')[1].split('.')[0])
		L = len(grid)

		# check if it exists
		for n in range(L):
			curlist = [grid[n][col] for col in ['folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
				'cav_theta','innerdustfile','outerdustfile','beta','L_sun','env_mass','env_power']]
			# fix a floating point issue
			for k in range(len(curlist)):
				if isinstance(curlist[k],float):
					curlist[k] = round(curlist[k],8)
					
			print "Curlist is:",curlist
			print "Mylist is:",mylist[1:]
			#print "Intersection:",set(mylist[1:]).intersection(set(curlist))
			#if len(set(mylist[1:]).intersection(set(curlist))) == len(set(curlist)):
			print curlist == [val for val in mylist[1:]]
			if curlist == [val for val in mylist[1:]]:
				add = 'no'
				print "Already found list in file, not adding"
				break
				
								
		# only add the new row when it's not already in there
		if add =='yes':
			print "Adding row"
			print mylist
			grid.add_row(mylist)

			# change all the names of the new lines
			grid['name'][-1] = grid['name'][-1]+"_"+str(Ngrid+1)
			print grid['name'][-1]
	
			i = -1
			model = m.YSOModelSim(name=grid['name'][i],folder=grid['folder'][i],L_sun=grid['L_sun'][i],M_sun=grid['M_sun'][i],env_mass=grid['env_mass'][i],\
				env_type='power',disk_rmax=grid['disk_rmax'][i],env_rmin=grid['env_rmin'][i],env_power=grid['env_power'][i],\
				disk_mass=grid['disk_mass'][i],cav=True,cav_theta=grid['cav_theta'][i],\
				amb_rmin=1,amb_rmax=grid['env_rmax'][i],env_rmax=grid['env_rmax'][i],T=grid['T'][i],innerdustfile=grid['innerdustfile'][i],\
				angles=angles,angles2=[0. for a in angles],env=True,disk=grid['disk'][i],\
				outerdustfile=grid['outerdustfile'][i],beta=grid['beta'][i])
			model.initModel()
			model.runModel()
			
			# save new grid after doing simulation
			f = open(filename,'wb')
			grid.write(folder[0]+name[0]+".grid.txt",format='ascii.fixed_width')
			pickle.dump(grid,f)
			f.close()

# this just applies to the first time we do a grid
else:	
	# this creates a grid in the form of a table
	allparams = ['name', 'folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
				'cav_theta','innerdustfile','outerdustfile','beta','L_sun','env_mass','env_power']
	types = ['<S30','<S30','f8','f8','f8','f8','<S30','f8','f8','f8','<S30','<S30','f8','f8','f8','f8']
	grid = Table(names=allparams,dtype=types)

	# only run things if the simulation has not yet been done
	for mylist in list_of_lists_flattened:
		grid.add_row(mylist)


	print "There will be ",len(grid),"models that will be run"
	for j in range(len(grid)):
		grid['name'][j] = grid['name'][j]+"_"+str(j)
	print grid
	grid.write(folder[0]+name[0]+".grid.txt",format='ascii.fixed_width')
	pickle.dump(grid,open(folder[0]+name[0]+".grid.dat",'wb'))

	for i in range(len(grid)):

		model = m.YSOModelSim(name=grid['name'][i],folder=grid['folder'][i],L_sun=grid['L_sun'][i],M_sun=grid['M_sun'][i],env_mass=grid['env_mass'][i],\
			env_type='power',disk_rmax=grid['disk_rmax'][i],env_rmin=grid['env_rmin'][i],env_power=grid['env_power'][i],\
			disk_mass=grid['disk_mass'][i],cav=True,cav_theta=grid['cav_theta'][i],\
			amb_rmin=1,amb_rmax=grid['env_rmax'][i],env_rmax=grid['env_rmax'][i],T=grid['T'][i],innerdustfile=grid['innerdustfile'][i],\
			angles=angles,angles2=[0. for a in angles],env=True,disk=grid['disk'][i],\
			outerdustfile=grid['outerdustfile'][i],beta=grid['beta'][i])
		model.initModel()
		model.runModel()
		#print grid[i]
		#model.plotModel(extinction=10,show=False,inc=4,sourcename=name[0],dist_pc=700) 	
		#model.plotSim(extinction=0,show=True,inc=inc,dist_pc=700)


#os.system("/cardini3/mrizzo/anaconda/bin/python run_chi_Squared.py")

