import models as m
import pickle
import time
from hyperion.util.convenience import OptThinRadius
'''
name = 'WL16-Ulrich-Env_Mass=0.1-Rc=2000-Lsun=5.6'
folder = 'WL16/'
WL16 = m.YSOModelSim(name=name,folder=folder,L_sun=5.6,env_mass=0.1,env_type='ulrich',rc=2000)
WL16.runModel()
WL16.plotModel(extinction=9.5,show=True)
WL16.modelDump()
'''
dustfilelist=['d03_5.5_3.0_A.hdf5']
dust_to_gas_ratio = 100
T=4400
env_rmin_list=[200]
env_rmax=28800
rc_list=[45]
l_sun_list=[30]#[25]
m_sun=5#2
env_mass_list=[1.]
mdisk=0.03
disk_rmax=100
h0 = 0.1
cav_theta=25.
disk_rmin=1
envtype='power'
env_power=-2
envmdot=1e-6
amb_rmin=1
amb_rmax=env_rmax

amb_dens = 6e-22
angles=[0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]
angles2=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
foldername='IRAS20050.3'
env=True	
cav=False
ext=10
inc=7
for env_rmin in env_rmin_list:
	for l_sun in l_sun_list:
		for env_mass in env_mass_list:
			for rc in rc_list:
				for dustfile in dustfilelist:
					if dustfile == 'd03_5.5_3.0_A.hdf5':
						mdisk /=dust_to_gas_ratio
						env_mass /=dust_to_gas_ratio
						envmdot /=dust_to_gas_ratio
					start = time.time()
					name = foldername+'-'+envtype+'-Env_Mass='+str(env_mass)+'-Rmin='+str(env_rmin)+'-Rmax='+str(env_rmax)+\
					'-Lsun='+str(l_sun)+'-Msun='+str(m_sun)+'-pow='+str(env_power)+'-h0='+str(h0)+'-mdisk='+str(mdisk)+ \
					'-rc='+str(rc)+'-cav='+str(cav_theta)+'-disk_rmin='+str(disk_rmin)+'-envmdot='+str(envmdot)+ \
					'-disk_rmax='+str(disk_rmax)+'-T='+str(T)+'-'+dustfile.split('.')[0]+'-angle='+str(angles[0])+ \
					'-ext='+str(ext)+'-amb_rmax='+str(amb_rmax)+'-draine'+'-OH5'
					if env==False:
						name += '-noenv'
					if cav==False:
						name += '-nocav'
					### 
					folder = foldername+'/'
					#WL16=m.modelLoad(folder,name)

					IRAS20050 = m.YSOModelSim(name=name,folder=folder,L_sun=l_sun,M_sun=m_sun,env_mass=env_mass,env_type=envtype,disk_rmax=disk_rmax,\
							env_rmin=env_rmin,rc=rc,env_power=env_power,disk_h_0=h0,disk_mass=mdisk,\
							cav=cav,cav_theta=cav_theta,disk_rmin=disk_rmin,mdot=envmdot,amb_rmin=amb_rmin,amb_rmax=amb_rmax,\
							env_rmax=env_rmax,T=T,innerdustfile=dustfile,amb_dens=amb_dens,angles=angles,angles2=angles2,env=env)
					IRAS20050.initModel()
					IRAS20050.runModel()
					IRAS20050.calcChi2(extinction=3,sourcename='IRAS20050.3',dist_pc=700)
					print 'Elapsed time: ',time.time()-start,'s'
					IRAS20050.plotModel(extinction=3,show=True,inc=inc,sourcename='IRAS20050.3',dist_pc=700) 	
					#time.sleep(2)
					#m.modelDump(WL16)


print 'Script is finished!'
