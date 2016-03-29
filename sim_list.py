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
T=5000
env_rmin_list=[300]
env_rmax=8000
rc_list=[700]
l_sun_list=[100]
env_mass_list=[0.05]
mdisk=0.00001
disk_rmax=450
h0 = 20
cav=65.
disk_rmin=1
envtype='ulrich'
env_power=0.0
envmdot=0.5e-6
amb_rmin=1
dustfilelist=['vsg.hdf5','OH5.hdf5']

foldername='WL16'

for env_rmin in env_rmin_list:
	for l_sun in l_sun_list:
		for env_mass in env_mass_list:
			for rc in rc_list:
				for dustfile in dustfilelist:
					start = time.time()
					name = folder+'-'+envtype+'-Env_Mass='+str(env_mass)+'-Rmin='+str(env_rmin)+'-Rmax='+str(env_rmax)+'-Lsun='+str(l_sun)+'-pow='+str(env_power)+'-h0='+str(h0)+'-mdisk='+str(mdisk)+'-rc='+str(rc)+'-cav='+str(cav)+'-disk_rmin='+str(disk_rmin)+'-envmdot='+str(envmdot)+'-disk_rmax='+str(disk_rmax)+'-noamb'+'-T='+str(T)+'-'+dustfile.split('.')[0]#+'-nodisk'
					### 
					folder = foldername+'/'
					#WL16=m.modelLoad(folder,name)

					WL16 = m.YSOModelSim(name=name,folder=folder,L_sun=l_sun,env_mass=env_mass,env_type=envtype,disk_rmax=disk_rmax,\
						env_rmin=env_rmin,rc=rc,env_power=env_power,disk_h_0=h0,disk_mass=mdisk,\
						cav=True,cav_theta=cav,disk_rmin=disk_rmin,mdot=envmdot,amb_rmin=amb_rmin,\
						env_rmax=env_rmax,T=T,dustfile=dustfile)
					WL16.runModel()
					print 'Elapsed time: ',time.time()-start,'s'
					WL16.plotModel(extinction=20,show=False) 	
					#time.sleep(2)
					#m.modelDump(WL16)


print 'Script is finished!'
