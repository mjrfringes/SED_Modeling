import models as m
import pickle
import subprocess as sp
import time
import os
#m.YSOModel(name='WL16-ResslerModel',folder='WL16/',env_mass=0.08)
#m.YSOModel(name='WL16-ResslerModel_no_envelope',folder='WL16/',env_mass=0.08,env=False,cav=False)

### Enoch values (bolometric luminosity with envelope)
#m.YSOModel(name='WL16-EnochModel',folder='WL16/',env_mass=0.08,L_sun=5.6)

#m.plotYSOModel(name='WL16-ResslerModel',folder='WL16/',extinction=31,show=False)
#m.plotYSOModel(name='WL16-ResslerModel_no_envelope',folder='WL16/',extinction=31,show=False)
#m.plotYSOModel(name='WL16-EnochModel',folder='WL16/',extinction=9.76,show=False)

###
'''
name = 'WL16-Test_Donut'
folder = 'WL16/'
WL16 = m.YSOModelSim(name=name,folder=folder,env_mass=0.08,env_type='ulrich',rc=400)
#WL16.runModel()
#WL16 = m.modelLoad(folder,name)
WL16.plotModel(extinction=9.5,show=True)
#WL16.modelDump()
print "Dumping file..."
try:
	os.system('rm %s' % (folder+name+'.mod'))
	pickle.dump(WL16,open(folder+name+'.mod','wb'))
except:
	print "Pickling failed :-("

WL16 = m.modelLoad(folder,name)
WL16.plotModel(extinction=9.5,show=True)

#name = 'WL16-Ulrich-Env_Mass=1-Rc=800-Lsun=5.6'
#folder = 'WL16/'
#WL16 = m.YSOModelSim(name=name,folder=folder,L_sun=5.6,env_mass=1,env_type='ulrich',rc=800)
#WL16.runModel()
#WL16.plotModel(extinction=9.5,show=True)
#WL16.modelDump()
#print "test3"
'''
while True:
	p=sp.Popen("/cardini3/mrizzo/anaconda/bin/python sim_list_WL6.py",stdout=sp.PIPE,shell=True)
	with p.stdout:
		for line in iter(p.stdout.readline,b''):
			print line,
	while p.poll() is None:
		time.sleep(1)
		print 'Process is not done (Ctrl-C to interrupt)'
	print 'Process is done - restarting the script (Ctrl-C to interrupt)'
#os.system('/cardini3/mrizzo/anaconda/bin/python sim_list.py')

