'''
This program runs hyperion on a grid or predefined models
'''
import itertools
from astropy.table import Table
import models as m
import pickle,os

### these are the parameters to probe:

disk = [1]
name = ['IRAS20050']
folder = ['IRAS20050/']
T_list=[4500.]
m_sun_list = [5.]
l_sun_list = [10.]
disk_mass_list = [0.003]
disk_rmax_list = [100.]
disk_rmin_list = [1.]
disk_h0_list = [0.1]


env_mass_list = [10.]
env_rmax_list = [10000.]
env_rmin_list = [10.,100.]

env_power_list = [-2.0]
env_type_list = ['power']
env_rc_list = [100.]
mdot_list = [1e-6]
cav_theta_list = [25.]

innerdust_list = ['OH5.hdf5']
outerdust_list = ['OH5.hdf5']

# this creates the list of lists
list_of_lists = [name, folder, l_sun_list, m_sun_list, env_mass_list, env_type_list,disk_rmax_list, env_rmin_list,env_rc_list,env_power_list, disk_h0_list, disk_mass_list , cav_theta_list, disk_rmin_list ,mdot_list, env_rmax_list, T_list,innerdust_list,outerdust_list,disk]

# flatten the list of lists
list_of_lists_flattened = list(itertools.product(*list_of_lists))

# if file exists, load existing grid file
filename = folder[0]+name[0]+'.grid.dat'
if os.path.exists(filename):
	grid = pickle.load(open(filename,'r'))
	# load up last line of grid
	Ngrid = int(grid['name'][-1].split('_')[1].split('.')[0])
	print Ngrid
	num = 0
	L=len(grid)
	for mylist in list_of_lists_flattened:	
		add = 'yes'
		for n in range(L):
			curlist = [grid[n][col] for col in ['folder','L_sun','M_sun','env_mass','env_type','disk_rmax','env_rmin','rc','env_power','disk_h0','disk_mass',\
				'cav_theta','disk_rmin','mdot','env_rmax','T','innerdustfile','outerdustfile','disk']]
			print "Curlist is:",set(mylist[1:])
			print "Mylist is:",set(curlist)
			print set(mylist[1:]).intersection(set(curlist))
			if len(set(mylist[1:]).intersection(set(curlist))) == len(curlist):
				add = 'no'
								
		# only add the new row when it's not already in there
		if add =='yes':
			print "Adding row"
			print mylist
			grid.add_row(mylist)
			num+=1
			print "Ngrid+num",Ngrid+num
			# change all the names of the new lines
			grid['name'][Ngrid+num] = grid['name'][Ngrid+num]+"_"+str(Ngrid+num)
			print grid['name'][Ngrid+num]
	
	# save new grid before doing simulation
	grid.write(folder[0]+name[0]+".grid.txt",format='ascii.fixed_width')
	pickle.dump(grid,open(folder[0]+name[0]+".grid.dat",'wb'))
	
	if num>0:
		for i in range(Ngrid+num,len(grid)):
			amb_dens = 1e-19
			if grid['innerdustfile'][i] == 'd03_5.5_3.0_A.hdf5':
				grid['disk_mass'][i] /= 100.
			if grid['outerdustfile'][i] == 'd03_5.5_3.0_A.hdf5':
				grid['env_mass'][i] /= 100.
				amb_dens /= 100.
	
	
			#if grid['grid['env_type'][i]'][i]==1:
			#	qdisk = True
			#else:
			#	qdisk=False

			model = m.YSOModelSim(name=grid['name'][i],folder=grid['folder'][i],L_sun=grid['L_sun'][i],M_sun=grid['M_sun'][i],env_mass=grid['env_mass'][i],\
				env_type=grid['env_type'][i],disk_rmax=grid['disk_rmax'][i],env_rmin=grid['env_rmin'][i],rc=grid['rc'][i],env_power=grid['env_power'][i],\
				disk_h_0=grid['disk_h0'][i],disk_mass=grid['disk_mass'][i],cav=True,cav_theta=grid['cav_theta'][i],disk_rmin=grid['disk_rmin'][i],\
				mdot=grid['mdot'][i],amb_rmin=1,amb_rmax=grid['env_rmax'][i],env_rmax=grid['env_rmax'][i],T=grid['T'][i],innerdustfile=grid['innerdustfile'][i],\
				amb_dens=amb_dens,angles=[0.,10.,20.,30.,40.,50.,60.,70.,80.,90.],angles2=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],env=True,disk=grid['disk'][i],\
				outerdustfile=grid['outerdustfile'][i])
			model.initModel()
			model.runModel()

# this just applies to the first time we do a grid
else:	
	# this creates a grid in the form of a table
	allparams = ['name', 'folder','L_sun','M_sun','env_mass','env_type','disk_rmax','env_rmin','rc','env_power','disk_h0','disk_mass',\
				'cav_theta','disk_rmin','mdot','env_rmax','T','innerdustfile','outerdustfile','disk']
	types = ['<S30','<S30','f8','f8','f8','<S30','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','<S30','<S30','i8']
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
		amb_dens = 1e-19
		if grid['innerdustfile'][i] == 'd03_5.5_3.0_A.hdf5':
			grid['disk_mass'][i] /= 100.
		if grid['outerdustfile'][i] == 'd03_5.5_3.0_A.hdf5':
			grid['env_mass'][i] /= 100.
			amb_dens /= 100.
	
	
		#if grid['grid['env_type'][i]'][i]==1:
		#	qdisk = True
		#else:
		#	qdisk=False

		model = m.YSOModelSim(name=grid['name'][i],folder=grid['folder'][i],L_sun=grid['L_sun'][i],M_sun=grid['M_sun'][i],env_mass=grid['env_mass'][i],\
			env_type=grid['env_type'][i],disk_rmax=grid['disk_rmax'][i],env_rmin=grid['env_rmin'][i],rc=grid['rc'][i],env_power=grid['env_power'][i],\
			disk_h_0=grid['disk_h0'][i],disk_mass=grid['disk_mass'][i],cav=True,cav_theta=grid['cav_theta'][i],disk_rmin=grid['disk_rmin'][i],\
			mdot=grid['mdot'][i],amb_rmin=1,amb_rmax=grid['env_rmax'][i],env_rmax=grid['env_rmax'][i],T=grid['T'][i],innerdustfile=grid['innerdustfile'][i],\
			amb_dens=amb_dens,angles=[0.,10.,20.,30.,40.,50.,60.,70.,80.,90.],angles2=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],env=True,disk=grid['disk'][i],\
			outerdustfile=grid['outerdustfile'][i])
		model.initModel()
		model.runModel()
		print grid[i]
		#model.plotModel(extinction=10,show=False,inc=4,sourcename=name[0],dist_pc=700) 	




