'''
plot multiple models
'''
import models as m
import pickle
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import seaborn as sns
from hyperion.model import ModelOutput
from hyperion.util.constants import rsun, lsun, au, msun, yr, c, pc, sigma
import photometry as phot

#model1 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_mdisk0.03_0.rtout"
#model2 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_0.rtout"
#model3 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_0.rtout"
#model4 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_mdisk0.01_0.rtout" # should be the same as model 2
#model5 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot10lsun_0.rtout"
#model6 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot10000lsun_0.rtout"
#model7 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot9lsun_0.rtout"
#model8 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot90lsun_0.rtout"
#model9 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot180lsun_0.rtout"
#model10 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_mdot2lsun1rmin0.1h0.rtout"
#model11 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_1rmin0.1h0_0.rtout"
#model12 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_1rmin0.1h0_coarse_.rtout"
#model13 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/flared_disk_0.1rmin0.01h0_coar.rtout"
#model14 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/flared_disk/alpha_disk_1rmin0.1h0coarse_0.rtout"
#model15 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/IRAS20050_0.rtout"

#model79 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/IRAS20050_79.rtout" #best IRA20050_2 fit with OH5 dust
#model711 ="/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/IRAS20050_711.rtout" #same as 79 but with draine dust
#model712 ="/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/IRAS20050_712.rtout" #same as 79 but with draine dust with a 100-1 gas-to-dust ratio
#model714 ="/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/IRAS20050_new/IRAS20050_714.rtout" #same as 79 but with kmh dust

#model18 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_18.rtout"

#model135 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_135.rtout"
#model60 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_60.rtout"
#model148 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_148.rtout"

model0 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_0.rtout"
model18 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_18.rtout"
model63 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_63.rtout"
model81 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_81.rtout"
model206 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_206.rtout"


model207 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_207.rtout"
model243 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_243.rtout" #disk mass x10


#IRAS2050.2
# best:
model27 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_27.rtout"
# mdisk x10
model244 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_244.rtout"
# rdisk /10
model245 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_245.rtout"
# renv /5
model246 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_246.rtout"
# renv *5
model247 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_247.rtout"
# renv_min = Tsub
model248 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_248.rtout"
# renv_min = Tsub & mdisk /100
model249 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_249.rtout"
# renv_min =10 and no disk
model250 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_250.rtout"
# same as 250 but with no ambient density
model296 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_296.rtout"

# envelope tests:
#normal:
model63= "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_63.rtout"
# envelope to 10000au instead of 5000 au
model307 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_307.rtout"
# envelope to 1000 au instead of 5000 au
model313 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_313.rtout"
# envelope to 100 au instead of 5000 au
model316 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_316.rtout"
# envelope to tiny mass
model317 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_317.rtout"
# envelope to tiny mass, disk /10
model318 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_318.rtout"
# envelope and disk both at 1e-4
model319 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_319.rtout"
# envelope and disk both at 1e-4; disk rmax to 10au, env rmin from 100au
model320 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_320.rtout"
# same as 320 but with no ambient density
model321 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_321.rtout"
# 1e-4 disk out to 100 au, envelope 3e-4
model322 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_322.rtout"
# 1e-4 disk out to 100 au, envelope 3e-3
model323 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_323.rtout"

# vary cavity opening angle
model346 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_346.rtout"
model347 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_347.rtout"
model348 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_348.rtout"
model349 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Grid/model_349.rtout"


#lsun= 3.846e+33
#2.103e+32 : alpha_disk_mdot9lsun1rmin
#1.928e+33 : alpha_disk_mdot9lsun0.1rmin 
## following with with 5 solar lum of source and 5 solar lum of accr (effectively, radius is 2 rsun)
#1.967e+34: alpha_disk_mdot5lsun0.05rmin0.005h0
#1.102e+34: alpha_disk_mdot5lsun0.1rmin0.01h0
## h0 doesn't change accretion luminosity
## rmin changes the accretion luminosity : 10 times smaller rmin means 10 times bigger lum
## following with with 8 solar lum of source and 2 solar lum of accr (effectively, radius is >2 rsun)
#1.511e+34: alpha_disk_mdot2lsun0.05rmin0.005h0
#4.694e+33: alpha_disk_mdot2lsun0.2rmin0.02h0
#1.037e+33: alpha_disk_mdot2lsun1rmin0.1h0
#7.087e+34:
#1160.75
model106 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_106.rtout" 
model107 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_107.rtout" 
model390 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_390.rtout" #107 w/ envelope size x2
model392 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_392.rtout" #107 w/ envelope size x10
model393 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_393.rtout" #exploring param
model395 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_395.rtout" #
model396 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_396.rtout" #
model397 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_397.rtout" # 
model398 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_398.rtout" # Lsun=400, Msun = 5 w/envelope of 5000
model399 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_399.rtout" # Lsun=100, Msun = 5 w/envelope of 5000
model400 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_400.rtout" # Lsun=100, Msun = 15 w/envelope of 5000
model401 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_401.rtout" # Lsun=100, Msun = 50 w/envelope of 5000
model402 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_402.rtout" # Lsun=200, Msun = 50 w/envelope of 5000
model403 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_403.rtout" # Lsun=200, Msun = 50 w/envelope of 5000; also put inner radius of envelope all the way down to OptRadius(1600): not really changing much
model404 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_404.rtout" # Lsun=200, Msun = 50 w/envelope of 5000; envelope power = -1.0 - not much change
model405 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_405.rtout" # Lsun=200, Msun = 50 w/envelope of 5000; envelope power = -1.0; no cavity - very interesting, one loses all of the short wavelength contrib!
model406 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_406.rtout" # Lsun=200, Msun = 50 w/envelope of 5000; envelope power = -1.0; cavity = 50deg --MISTAKE--
model407 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_407.rtout" # Lsun=200, Msun = 50 w/envelope of 5000; envelope power = -1.0; cavity = 50deg 
model408 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_408.rtout" # Lsun=120, Msun = 50 w/envelope of 5000; envelope power = -1.0; cavity = 50deg 
model409 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_409.rtout" # Lsun=120, Msun = 10 w/envelope of 5000; envelope power = -1.0; cavity = 50deg, disk_rmanx=200, disk_mass=0.2, env=10
model410 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_410.rtout" # Lsun=120, Msun = 10 w/envelope of 5000; envelope power = -2.0; cavity = 50deg, disk_rmanx=200, disk_mass=0.2, env=10
model411 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_411.rtout" # Lsun=120, Msun = 10 w/envelope of 5000; envelope power = -2.0; cavity = 50deg, disk_rmanx=200, disk_mass=0.2, env=10, OH5 outside
model412 = "/cardini3/mrizzo/2012SOFIA/SED_Models/hyperion/Gridfinal/model_412.rtout" # Lsun=120, Msun = 15 w/envelope of 5000; envelope power = -2.0; cavity = 50deg, disk_rmanx=200, disk_mass=0.2, env=10, OH5 everywhere, 10 times more photons

#model_list = [model63,model313,model322]
model_list = [model410,model411,model412]
sns.set(context='paper',style='ticks')
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.labelsize']=16
mpl.rcParams['legend.fontsize']=14
mpl.rcParams['font.size']=14

def distance(target):
	if "IRAS20050" in target:
		return 700.*pc
	elif "NGC2071" in target:
		return 420.*pc

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))
sources = sourcetable.group_by('SOFIA_name')

markers = ['v','p','D','^','h','o','*','>','<']
TwoMASS = phot.plots(['j','h','ks'],[1.3,1.6,2.2],sns.xkcd_rgb['pale red'],markers[0],'2MASS')
Spitzer = phot.plots(['i1','i2','i3','i4','m1','m2'],[3.6,4.5,5.8,8.,24,70],sns.xkcd_rgb['medium green'],markers[1],'Spitzer')
SOFIA = phot.plots(['F11','F19','F31','F37'],[11.1,19.7,31.5,37.1],sns.xkcd_rgb['denim blue'],markers[3],'SOFIA')
HERSCHEL = phot.plots(['H70','H160','H250','H350'],[70.,160.,250.,350.],sns.xkcd_rgb['purple'],markers[4],'HERSCHEL')
VANKEMPEN = phot.plots(['SMA1300'],[1300],sns.xkcd_rgb['denim blue'],markers[5],'VANKEMPEN')
fluxlist = [TwoMASS,Spitzer,HERSCHEL,SOFIA,VANKEMPEN]
fluxnames = [p.bands for p in fluxlist]


fig = plt.figure(figsize=(15,7.5))
ax = fig.add_subplot(111)
colors = [sns.xkcd_rgb['denim blue'],sns.xkcd_rgb['pale red'],sns.xkcd_rgb['medium green']] #sns.color_palette('husl',9)
target="NGC2071.2"
for i in range(len(model_list)):
		model = model_list[i]
		mo = ModelOutput(model)
		sed = mo.get_sed(aperture=-1, inclination='all', distance=distance(target))
		ax.plot(sed.wav,sed.val.transpose(),color=colors[i],alpha=0.8)
		
for key,source in zip(sources.groups.keys,sources.groups):
	if target == source['SOFIA_name'][0]:	
		for fl in fluxlist:
			#`print fl.bands
			phot.plotData(ax,source,fl,0.8,msize=12)


ax.set_xscale('log')
ax.set_yscale('log')

plt.show()

