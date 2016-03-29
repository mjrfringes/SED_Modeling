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
T=3000
env_rmin_list=[0.25]
env_rmax=1300
rc_list=[76]
l_sun_list=[10]#[25]
m_sun=0.2#2
env_mass_list=[0.05]
mdisk=0.00002
disk_rmax=70
h0 = 20
cav_theta=42.
disk_rmin=0.1
envtype='ulrich'
env_power=0.0
envmdot=0.4e-6
amb_rmin=1
amb_rmax=env_rmax
dustfilelist=['OH5.hdf5']
amb_dens = 2e-21
angles=[0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]
angles2=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
foldername='SR21'
env=True	
cav=True
ext=3
inc=0
for env_rmin in env_rmin_list:
	for l_sun in l_sun_list:
		for env_mass in env_mass_list:
			for rc in rc_list:
				for dustfile in dustfilelist:
					start = time.time()
					name = foldername+'-'+envtype+'-Env_Mass='+str(env_mass)+'-Rmin='+str(env_rmin)+'-Rmax='+str(env_rmax)+\
					'-Lsun='+str(l_sun)+'-Msun='+str(m_sun)+'-pow='+str(env_power)+'-h0='+str(h0)+'-mdisk='+str(mdisk)+ \
					'-rc='+str(rc)+'-cav='+str(cav_theta)+'-disk_rmin='+str(disk_rmin)+'-envmdot='+str(envmdot)+ \
					'-disk_rmax='+str(disk_rmax)+'-T='+str(T)+'-'+dustfile.split('.')[0]+'-angle='+str(angles[0])+ \
					'-ext='+str(ext)+'-amb_rmax='+str(amb_rmax)+'-kmh' #+'-nodisk'
					### 
					folder = foldername+'/'
					#WL16=m.modelLoad(folder,name)

					SR21 = m.YSOModelSim(name=name,folder=folder,L_sun=l_sun,M_sun=m_sun,env_mass=env_mass,env_type=envtype,disk_rmax=disk_rmax,\
							env_rmin=env_rmin,rc=rc,env_power=env_power,disk_h_0=h0,disk_mass=mdisk,\
							cav=cav,cav_theta=cav_theta,disk_rmin=disk_rmin,mdot=envmdot,amb_rmin=amb_rmin,amb_rmax=amb_rmax,\
							env_rmax=env_rmax,T=T,dustfile=dustfile,amb_dens=amb_dens,angles=angles,angles2=angles2,env=env)
					SR21.runModel()
					print 'Elapsed time: ',time.time()-start,'s'
					SR21.plotModel(extinction=ext,show=True,inc=0,sourcename='Oph.1') 	
					#time.sleep(2)
					#m.modelDump(WL16)


print 'Script is finished!'
