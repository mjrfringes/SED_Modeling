### this is a python script that runs simples models to determine the spatial extension of WL16


### dependencies ###
import numpy as np
import matplotlib.pyplot as plt
from hyperion.dust import HenyeyGreensteinDust,IsotropicDust,SphericalDust
from hyperion.model import AnalyticalYSOModel,ModelOutput
from hyperion.util.constants import rsun, lsun, au, msun, yr, c, pc, sigma
from hyperion.util.convenience import OptThinRadius
from hyperion.grid import SphericalPolarGrid
import astropy.constants as const
from astropy.io import fits
from scipy.interpolate import interp1d
import pickle
import os
import subprocess as sp
import time

#import sys
#sys.path.append('/cardini3/mrizzo/SOFIA2012/Reduction')
import photometry as p
import seaborn as sns


colors = sns.color_palette('hls',9) ### this picks the color palette

class YSOModelSim(object):
	
	def __init__(self,name,folder,T=9000,M_sun=5.6,L_sun=250,disk_mass=0.01,disk_rmax=100, 
		env=True,env_type='power',rc=400,mdot=1e-8,env_mass=0.1,env_rmin=30,env_rmax=5000,cav=True,cav_r0=500,cav_rho_0=8e-24,cav_theta=25,env_power=-1.5,
		Npix=149,angles=[20.,45.,60.,80],angles2=[60.,60.,60.,60.], amb_dens=8e-24, disk="Flared",disk_rmin=1., amb_rmin=1., amb_rmax=1000., innerdustfile='OH5.hdf5',
		outerdustfile='d03_5.5_3.0_A.hdf5',beta=1.1):
		self.name=name
		self.folder=folder
		self.T=T
		self.M_sun=M_sun*msun
		self.L_sun=L_sun*lsun
		self.disk_mass=disk_mass*msun
		self.disk_rmax=disk_rmax*au
		self.disk_rmin=disk_rmin*au
		self.disk_h_0 = OptThinRadius(1600)
		self.env=env
		self.disk=disk
		self.env_type=env_type
		self.env_mass=env_mass*msun
		self.env_rmin=env_rmin*au
		self.env_rmax=env_rmax*au
		self.mdot=mdot #*msun/yr*self.M_sun # disk accretion rate
		self.rc=rc*au
		self.cav=cav
		self.cav_rho_0=cav_rho_0
		self.cav_r0=cav_r0*au
		self.cav_theta=cav_theta
		self.Npix=Npix
		self.angles=angles
		self.angles2=angles2
		self.amb_dens=amb_dens
		self.amb_rmin=amb_rmin
		self.amb_rmax=amb_rmax*au
		self.env_power=env_power
		self.dustfile=innerdustfile
		self.dustfile_out=outerdustfile
		self.limval = max(self.env_rmax,1000*au)
		self.beta = beta

	def modelDump(self):
		sp.call('rm %s.mod ' % (self.folder+self.name),shell=True)
		pickle.dump(self,open(self.folder+self.name+'.mod','wb'))
		time.sleep(2)

	def modelPrint(self):
		#string= self.folder+ self.name+'\n'
		string="T="+str(self.T)+"K"+'\n'
		string+= "M="+str(self.M_sun/msun)+'Msun'+'\n'
		string+= "L="+str(self.L_sun/lsun)+'Lsun'+'\n'
		string+= "Disk="+str(self.disk)+'\n'
		string+= "Disk_mass="+str(self.disk_mass/msun)+'Msun'+'\n'
		string+= "Disk_rmax="+str(self.disk_rmax/au)+'AU'+'\n'
		string+= "Disk_rmin="+str(self.disk_rmin/au)+'AU'+'\n'
		string+= "env="+str(self.env)+'\n'
		string+= "env_type="+self.env_type+'\n'
		string+= "env_mass="+str(self.env_mass/msun)+'Msun'+'\n'
		string+= "env_rmax="+str(self.env_rmax/au)+'AU'+'\n'
		string+= "env_rmin="+str(self.env_rmin/au)+'AU'+'\n'
		if self.env_type == 'ulrich' and self.env==True:
			string+= "mass_ulrich="+str((8.*np.pi*self.env_rho_0*self.rc**3*pow(self.env_rmax/self.rc,1.5)/(3.*np.sqrt(2)))/msun)+'Msun'+'\n'
		string+= "mdot="+str(self.mdot)+'Msun/yr'+'\n' # (only if env_type="Ulrich")
		string+= "rc="+str(self.rc/au)+'AU'+'\n' # (only if env_type="Ulrich")
		string+= "cav="+str(self.cav)+'\n'
		string+= "cav_theta="+str(self.cav_theta)+'\n'
		string+= "cav_r0="+str(self.cav_r0/au)+'\n'
		string+= "env_power="+str(self.env_power)+'\n'
		string+= "disk_h_0="+str(self.disk_h_0)+'\n'
		string+= "dustfile="+self.dustfile+'\n'
		string+= "dustfile_out="+self.dustfile_out+'\n'
		string+= "amb_dens="+str(self.amb_dens)+'\n'
		string+= "amb_rmin="+str(self.amb_rmin)+'\n'
		string+= "amb_rmax="+str(self.amb_rmax/au)+'\n'
		string+= "angles="+str(self.angles)+'\n'
		print string
		return string

	def dust_gen(self,dustfile,dustfile_out='d03_5.5_3.0_A.hdf5'):
		### first, we need to load Tracy's dust files and manipulate them to feed to Hyperion
		### wavelength (microns),Cext,Csca,Kappa,g,pmax,theta (ignored)
		### albedo = Csca/Cext
		### opacity kappa is in cm^2/gm, dust_gas extinction opactity (absorption+scattering) - assumes gas-to=dust raio of 100
		### see Whitney et al. 2003a
		
#		tracy_dust = np.loadtxt('Tracy_models/OH5.par')

#		### format for dust: d = HenyeyGreensteinDust(nu, albedo, chi, g, p_lin_max)
#		nu = const.c.value/ (tracy_dust[:,0]*1e-6)
#		albedo = tracy_dust[:,2]/tracy_dust[:,1]
#		chi = tracy_dust[:,3]
#		g = tracy_dust[:,4]
#		p_lin_max = tracy_dust[:,5]

#		### flip the table to have an increasing frequency
#		nu = nu[::-1]
#		albedo = albedo[::-1]
#		chi=chi[::-1]
#		g=g[::-1]
#		p_lin_max=p_lin_max[::-1]

#		### create the dust model
#		d = HenyeyGreensteinDust(nu, albedo, chi, g, p_lin_max)
#		d.optical_properties.extrapolate_wav(0.001,1.e7)
#		d.plot('OH5.png')
#		d.write('OH5.hdf5')
		
		self.d = SphericalDust()
		self.d.read(dustfile)
		self.d.plot(str(dustfile.split(',')[:-1])+'.png')
		self.d_out = SphericalDust()
		self.d_out.read(dustfile_out)
		#self.d_out.read(dustfile)
		self.d_out.plot(str(dustfile_out.split(',')[:-1])+'.png')

	def initModel(self):
		### Use Tracy parameter file to set up the model 
		self.dust_gen(self.dustfile,self.dustfile_out)
		mi = AnalyticalYSOModel()

		mi.star.temperature = self.T
		mi.star.mass = self.M_sun
		mi.star.luminosity = self.L_sun
		mi.star.radius=np.sqrt(mi.star.luminosity/(4.0*np.pi*sigma*mi.star.temperature**4))
		#m.star.luminosity = 4.0*np.pi*m.star.radius**2*sigma*m.star.temperature**4
		print mi.star.luminosity/lsun
		self.luminosity=mi.star.luminosity/lsun

		if self.disk=="Flared":
			print "Adding flared disk"
			disk = mi.add_flared_disk()
			disk.dust=self.d
			if self.dustfile == 'd03_5.5_3.0_A.hdf5':
				disk.mass=self.disk_mass/100.
			else: disk.mass=self.disk_mass
			disk.rmin=OptThinRadius(1600) #self.disk_rmin
			print "disk.rmin = ",disk.rmin,disk.rmin/au
			disk.rmax=self.disk_rmax
			disk.r_0 = self.disk_rmin
			disk.h_0 = disk.r_0/10. #self.disk_h_0*au
			disk.beta=self.beta
			disk.p = -1.
		elif self.disk=="Alpha":
			print "Adding alpha disk"
			disk = mi.add_alpha_disk()
			disk.dust=self.d
			if self.dustfile == 'd03_5.5_3.0_A.hdf5':
				disk.mass=self.disk_mass/100.
			else: disk.mass=self.disk_mass
			disk.rmin=OptThinRadius(1600)
			disk.rmax=self.disk_rmax
			disk.r_0 = self.disk_rmin
			disk.h_0 = disk.r_0/10. #self.disk_h_0*au
			disk.beta=1.1
			disk.p = -1
			disk.mdot=self.mdot
			disk.star = mi.star
			
		#print 'Disk density:',disk.rho_0

		
		if self.env==True and self.env_type=='power':
			envelope=mi.add_power_law_envelope()
			envelope.dust=self.d_out
			envelope.r_0=self.env_rmin
			#envelope.r_0 = OptThinRadius(1600)
			if self.dustfile_out == 'd03_5.5_3.0_A.hdf5':
				envelope.mass=self.env_mass/100.
			else: envelope.mass=self.env_mass
			envelope.rmin=self.env_rmin
			envelope.rmax=self.env_rmax
			envelope.power=self.env_power
			#print 'Envelope rho:',envelope.rho_0
		elif self.env==True and self.env_type=='ulrich':
			envelope=mi.add_ulrich_envelope()
			envelope.dust=self.d_out
			envelope.mdot=1e-6*msun/yr # has little impact on the fluxes, so fixed
			envelope.rc=self.rc
			envelope.rmin=self.env_rmin
			envelope.rmax=self.env_rmax
		if self.env==True:
			self.env_rho_0 = envelope.rho_0
			print 'Envelope rho:',envelope.rho_0

		#print "Rho_0 = ",envelope.rho_0
		if self.cav==True:
			cavity=envelope.add_bipolar_cavity()
			cavity.dust=self.d_out
			cavity.power=1.5
			cavity.cap_to_envelope_density=True ### prevents the cavity density to go above the envelope's density
			cavity.r_0=self.cav_r0
			cavity.theta_0=self.cav_theta
			cavity.rho_0=self.cav_rho_0 #in g/cm^3
			cavity.rho_exp=0.0
			
		
#		if self.env==True:
#			ambient=mi.add_ambient_medium(subtract=[envelope,disk])
#		if self.dustfile_out == 'd03_5.5_3.0_A.hdf5':
#			ambient.rho=self.amb_dens/100.
#		else: ambient.rho=self.amb_dens
#		ambient.rmin=OptThinRadius(1600.)
#		ambient.rmax=self.env_rmax
#		ambient.dust=self.d_out
		

		'''*** Grid parameters ***'''
		mi.set_spherical_polar_grid_auto(199,49,1)

		# Specify that the specific energy and density are needed
		mi.conf.output.output_specific_energy = 'last'
		mi.conf.output.output_density = 'last'


		'''**** Output Data ****'''
		image = mi.add_peeled_images(sed=True,image=False)
		image.set_wavelength_range(150,1,3000)
		#image.set_image_size(self.Npix,self.Npix)
		#image.set_image_limits(-self.limval,self.limval,-self.limval,self.limval)
		image.set_aperture_range(1,100000.*au,100000.*au)
		image.set_viewing_angles(self.angles,self.angles2)
		#image.set_track_origin('detailed')
		image.set_uncertainties(True)

		''' Use the modified random walk
		*** Advanced ***'
		YES = DIFFUSION  = Whether to use the diffusion
		'''
		if self.env==True:
			#mi.set_pda(True)
			mi.set_mrw(True)
		else:
			mi.set_pda(False)
			mi.set_mrw(False)

		# Use raytracing to improve s/n of thermal/source emission
		mi.set_raytracing(True)


		'''**** Preliminaries ****'''
		mi.set_n_initial_iterations(5)
		mi.set_n_photons(initial=1e6,imaging=1e6,raytracing_sources=1e5,raytracing_dust=1e6)
		mi.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)
		self.m = mi

	def runModel(self):
		self.initModel()
		self.m.write(self.folder+self.name+'.rtin')
		self.m.run(self.folder+self.name+'.rtout', mpi=True,n_processes=6)

	def plotData(self,ax,sourcename):
		if sourcename != 'None':
			folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
			sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

			markers = ['v','p','D','^','h','o','*','x','d','<']
			TwoMASS = ['j','h','ks']
			uTwoMASS = ["e_"+col for col in TwoMASS]
			wlTwoMASS = [1.3,1.6,2.2]
			colTwoMASS = colors[0]
			markerTwoMASS = markers[0]
			labelTwoMASS = '2MASS'
			Spitzer = ['i1','i2','i3','i4','m1','m2']
			uSpitzer = ["e_"+col for col in Spitzer]
			wlSpitzer = [3.6,4.5,5.8,8.,24,70]
			colSpitzer = colors[1]
			markerSpitzer = markers[1]
			labelSpitzer = 'Spitzer'
			WISE = ['w1','w2','w3','w4']
			uWISE = ["e_"+col for col in WISE]
			wlWISE = [3.4,4.6,12,22]
			colWISE = colors[2]
			labelWISE = 'WISE'
			markerWISE = markers[2]
			SOFIA = ['F11','F19','F31','F37']
			uSOFIA = ["e_"+col for col in SOFIA]
			wlSOFIA = [11.1,19.7,31.5,37.1]
			colSOFIA = colors[3]
			markerSOFIA = markers[3]
			labelSOFIA = 'SOFIA'
			IRAS = ['Fnu_12','Fnu_25','Fnu_60','Fnu_100']
			uIRAS = ["e_"+col for col in IRAS]
			wlIRAS = [12,25,60,100]
			colIRAS = colors[4]
			markerIRAS = markers[4]
			labelIRAS = 'IRAS'
			AKARI = ['S65','S90','S140','S160']
			uAKARI = ["e_"+col for col in AKARI]
			wlAKARI = [65,90,140,160]
			colAKARI = colors[5]
			markerAKARI = markers[5]
			labelAKARI = 'AKARI'
			ENOCH = ['Fp']
			uENOCH = ["e_"+col for col in ENOCH]
			wlENOCH = [1300]
			colENOCH = colors[6]
			markerENOCH = markers[6]
			labelENOCH = 'ENOCH'
			HERSCHEL = ['H70','H160','H250','H350','H500']
			uHERSCHEL = ["e_"+col for col in HERSCHEL]
			wlHERSCHEL = [70,160,250,350,500]
			colHERSCHEL = colors[7]
			markerHERSCHEL = markers[7]
			labelHERSCHEL = 'HERSCHEL'
			SCUBA = ['S450','S850','S1300']
			uSCUBA = ["e_"+col for col in SCUBA]
			wlSCUBA = [450,850,1300]
			colSCUBA = colors[8]
			markerSCUBA = markers[8]
			labelSCUBA = 'SCUBA'
			alpha=1
			sources = sourcetable.group_by('SOFIA_name')
			for key,sourcetable in zip(sources.groups.keys,sources.groups):
				if sourcename == sourcetable['SOFIA_name'][0]:	
					#print sourcetable['SOFIA_name'][0]
					p.plotData(ax,sourcetable,markerTwoMASS,TwoMASS,uTwoMASS,wlTwoMASS,colTwoMASS,labelTwoMASS,alpha)
					p.plotData(ax,sourcetable,markerSpitzer,Spitzer,uSpitzer,wlSpitzer,colSpitzer,labelSpitzer,alpha)
					p.plotData(ax,sourcetable,markerWISE,WISE,uWISE,wlWISE,colWISE,labelWISE,alpha)
					p.plotData(ax,sourcetable,markerSOFIA,SOFIA,uSOFIA,wlSOFIA,colSOFIA,labelSOFIA,alpha)
					p.plotData(ax,sourcetable,markerIRAS,IRAS,uIRAS,wlIRAS,colIRAS,labelIRAS,alpha)
					p.plotData(ax,sourcetable,markerAKARI,AKARI,uAKARI,wlAKARI,colAKARI,labelAKARI,alpha)
					p.plotData(ax,sourcetable,markerENOCH,ENOCH,uENOCH,wlENOCH,colENOCH,labelENOCH,alpha)
					p.plotData(ax,sourcetable,markerHERSCHEL,HERSCHEL,uHERSCHEL,wlHERSCHEL,colHERSCHEL,labelHERSCHEL,alpha)
					p.plotData(ax,sourcetable,markerSCUBA,SCUBA,uSCUBA,wlSCUBA,colSCUBA,labelSCUBA,alpha)

	def calcChi2(self,dist_pc=140,extinction=0, sourcename='Oph.1'):
		self.dist=dist_pc*pc
		self.extinction=extinction
		chi = np.loadtxt('kmh94_3.1_full.chi')
		wav = np.loadtxt('kmh94_3.1_full.wav')
		Chi = interp1d(wav,chi,kind='linear')
		modelname = self.folder+self.name
		self.mo = ModelOutput(modelname+'.rtout')
		
		# get the sed of all inclination
		sed = self.mo.get_sed(aperture=-1, inclination='all', distance=self.dist,units='Jy')
				
		# calculate the optical depth at all wavelengths
		tau = self.extinction*Chi(sed.wav)/Chi(0.550)/1.086
		
		# calculate extinction values
		ext = np.array([np.exp(-tau) for i in range(sed.val.shape[0])])
		
		# apply extinction to model
		extinct_values = np.log10(sed.val.transpose()*ext.T)
		
		# data points and errors
		folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
		sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))
		TwoMASS = ['j','h','ks']
		uTwoMASS = ["e_"+col for col in TwoMASS]
		wlTwoMASS = [1.3,1.6,2.2]
		labelTwoMASS = '2MASS'
		Spitzer = ['i1','i2','i3','i4']
		uSpitzer = ["e_"+col for col in Spitzer]
		wlSpitzer = [3.6,4.5,5.8,8.]
		labelSpitzer = 'Spitzer'
		SOFIA = ['F11','F19','F31','F37']
		uSOFIA = ["e_"+col for col in SOFIA]
		wlSOFIA = [11.1,19.7,31.5,37.1]
		labelSOFIA = 'SOFIA'
		sources = sourcetable.group_by('SOFIA_name')
		for key,source in zip(sources.groups.keys,sources.groups):
			if sourcename == source['SOFIA_name'][0]:	
				datapoints = source[TwoMASS+Spitzer+SOFIA]
				dataerrors = source[uTwoMASS+uSpitzer+uSOFIA]
				print p.nptable(datapoints),p.nptable(dataerrors)
				
				# calculate log10 of quantities required for chi squared
				logFnu = np.log10(p.nptable(datapoints))-0.5*(1./np.log(10.))*p.nptable(dataerrors)**2/p.nptable(datapoints)**2
				varlogFnu = (1./np.log(10)/p.nptable(datapoints))**2*p.nptable(dataerrors)**2
				print extinct_values,extinct_values.shape
				
				# for each inclination, calculate chi squared; need to interpolate to get model at required wavelengths
				Ninc = extinct_values.shape[1]
				chi2 = np.zeros(Ninc)
				wl=wlTwoMASS+wlSpitzer+wlSOFIA
				N = len(wl)
				for j in range(Ninc):
					interp_func = interp1d(sed.wav,extinct_values[:,j],kind='linear')
					interp_vals = interp_func(wl)
					chi2[j] = 1./N * np.sum((logFnu - interp_vals)**2/varlogFnu)
					
				print chi2

	def plotModel(self,dist_pc=140,inc=3,extinction=0,show=False,sourcename='Oph.1'):
		self.dist=dist_pc*pc
		self.inc=inc
		self.extinction=extinction
		modelname = self.folder+self.name
		self.mo = ModelOutput(modelname+'.rtout')

		#tracy_dust = np.loadtxt('Tracy_models/OH5.par')
		chi = np.loadtxt('kmh94_3.1_full.chi')
		wav = np.loadtxt('kmh94_3.1_full.wav')
		Chi = interp1d(wav,chi,kind='linear')



		fig = plt.figure(figsize=(20,14))
		ax=fig.add_subplot(2,3,1)
		sed = self.mo.get_sed(aperture=-1, inclination='all', distance=self.dist)
		#print tracy_dust[11,1],Cext(sed.wav[-1]),Cext(sed.wav[-1])/tracy_dust[11,1]
		tau = self.extinction*Chi(sed.wav)/Chi(0.550)/1.086
		#print Cext(sed.wav)/tracy_dust[11,1]
		ext = np.array([np.exp(-tau) for i in range(sed.val.shape[0])])
		#print tau,np.exp(-tau)
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='black')
		ax.set_title(modelname+'_seds, Av='+str(self.extinction))
		ax.set_xlim(sed.wav.min(), 1300)
		ax.set_ylim(1e-13, 1e-7)
		ax.set_xlabel(r'$\lambda$ [$\mu$m]')
		ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
		self.plotData(ax,sourcename)
		ax.set_xscale('log')
		ax.set_yscale('log')

		#ax.set_ylabel(r'$F_{Jy}$ [Jy]')
		#plt.legend(loc=4)

		ax=fig.add_subplot(2,3,2)
		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist)
		ext=np.exp(-tau)
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, lw=3,color='black',label='source_total')
		ax.set_xlim(sed.wav.min(), 1300)
		ax.set_ylim(1e-13, 1e-7)  ### for lamFlam
		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='source_emit')
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='blue',label='source_emit')
		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='source_scat')
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='teal',label='source_scat')
		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='dust_emit')
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='red',label='dust_emit')
		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='dust_scat')
		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='orange',label='dust_scat')
		self.plotData(ax,sourcename)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_title('seds_inc=inc')
		ax.set_xlabel(r'$\lambda$ [$\mu$m]')
		ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
		#ax.set_ylabel(r'$F_{Jy}$ [Jy]')
		leg = ax.legend(loc=4,fontsize='small')
		#leg = plt.gca().get_legend()
		#plt.setp(leg.get_text(),fontsize='small')
		# Extract the quantities
		g = self.mo.get_quantities()
		
		# Get the wall positions for r and theta
		rw, tw = g.r_wall / au, g.t_wall

		# Make a 2-d grid of the wall positions (used by pcolormesh)
		R, T = np.meshgrid(rw, tw)

		# Calculate the position of the cell walls in cartesian coordinates
		X, Z = R * np.sin(T), R * np.cos(T)

		# Make a plot in (x, z) space for different zooms
		from matplotlib.colors import LogNorm,PowerNorm
		# Make a plot in (r, theta) space
		ax = fig.add_subplot(2, 3, 3)
		if g.shape[-1]==2:
			c = ax.pcolormesh(X, Z, g['temperature'][0].array[0, :, :]+g['temperature'][1].array[0, :, :],norm=PowerNorm(gamma=0.5,vmin=1,vmax=500))
		else :
			c = ax.pcolormesh(X, Z, g['temperature'][0].array[0, :, :],norm=PowerNorm(gamma=0.5,vmin=1,vmax=500))
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlim(X.min(), X.max()/5.)
		ax.set_ylim(Z.min()/10., Z.max()/10.)
		ax.set_xlabel('x (au)')
		ax.set_ylabel('z (au)')
		#ax.set_yticks([np.pi, np.pi * 0.75, np.pi * 0.5, np.pi * 0.25, 0.])
		#ax.set_yticklabels([r'$\pi$', r'$3\pi/4$', r'$\pi/2$', r'$\pi/4$', r'$0$'])
		cb = fig.colorbar(c)
		ax.set_title('Temperature structure')
		cb.set_label('Temperature (K)')
		#fig.savefig(modelname+'_temperature_spherical_rt.png', bbox_inches='tight')


		ax = fig.add_subplot(2, 3, 4)
		if g.shape[-1]==2:
			c = ax.pcolormesh(X, Z, g['density'][0].array[0, :, :]+g['density'][1].array[0, :, :],norm=LogNorm(vmin=1e-22,vmax=g['density'][0].array[0, :, :].max()))
		else :
			c = ax.pcolormesh(X, Z, g['density'][0].array[0, :, :],norm=LogNorm(vmin=1e-22,vmax=g['density'][0].array[0, :, :].max()))
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlim(X.min(), X.max()/5.)
		ax.set_ylim(Z.min()/10., Z.max()/10.)
		ax.set_xlabel('x (au)')
		ax.set_ylabel('z (au)')
		ax.set_title('Density structure')
		cb = fig.colorbar(c)
		cb.set_label('Density (g/cm2)')

		### plot the convolved image with the 37 micron filter (manually set to slice 18 of the cube - this would change with wavelength coverage)
		ax = fig.add_subplot(2, 3, 5)
		self.image = self.mo.get_image(inclination=inc,distance=self.dist,units='Jy')
		fits.writeto(modelname+'_inc_'+str(inc)+'.fits',self.image.val.swapaxes(0,2).swapaxes(1,2),clobber=True)

		### need to convolve the image with a Gaussian PSF
		pix = 2.*self.limval/au/self.Npix # in AU/pix
		pix_asec = pix/(self.dist/pc) # in asec/pix
		airy_asec = 3.5 #asec
		airy_pix = airy_asec/pix_asec # in pix
		gauss_pix = airy_pix/2.35 # in Gaussian std 
		print "Gaussian std: ",gauss_pix

		from scipy.ndimage.filters import gaussian_filter as gauss
		#print [(i,sed.wav[i]) for i in range(len(sed.wav))]

		img37 = self.image.val[:,:,18]
		convol = gauss(img37,gauss_pix,mode='constant',cval=0.0)
		Nc = self.Npix/2
		hw = min(int(20./pix_asec),Nc) #(max is Nc)
		#ax.imshow(img37,norm=LogNorm(vmin=1e-20,vmax=img37.max()))
		#ax.imshow(img37,interpolation='nearest')
		#ax.imshow(convol,norm=LogNorm(vmin=1e-20,vmax=img37.max()))
		#ax.imshow(convol,interpolation='nearest',norm=LogNorm(vmin=1e-20,vmax=img37.max()))
		ax.imshow(convol[Nc-hw:Nc+hw,Nc-hw:Nc+hw],interpolation='nearest',origin='lower',cmap=plt.get_cmap('gray'))
		airy_disk = plt.Circle((airy_pix*1.3,airy_pix*1.3),airy_pix,color=colors[3])		
		ax.add_artist(airy_disk)
		ax.text(airy_pix*3,airy_pix*1.3/2.0,'SOFIA 37um Airy disk',color=colors[3])
		ax.set_title('Convolved image')
		fits.writeto(modelname+'_inc_'+str(inc)+'_convol37.fits',convol,clobber=True)

		### draw a cross-section of the image to show the spatial extension in linear scale, to compare with what we observe in the model.
		ax = fig.add_subplot(2, 3, 6)
		ax.plot(range(Nc-hw,Nc+hw),convol[Nc-hw:Nc+hw,Nc-1],label='cross-section 1')
		ax.plot(range(Nc-hw,Nc+hw),convol[Nc-1,Nc-hw:Nc+hw],label='cross-section 2')
		maxconvol = convol[Nc-hw:Nc+hw,Nc-1].max()
		gauss = np.exp( -(np.array(range(-hw,hw))**2 / (2. * gauss_pix**2)))
		gauss/= gauss.max()
		gauss*=maxconvol
		ax.plot(range(Nc-hw,Nc+hw),gauss,label='SOFIA beam')
		leg = ax.legend(loc=2,fontsize='small')
		#leg = plt.gca().get_legend()
		#plt.setp(leg.get_text(),fontsize='small')
		ax.set_title('Cross section at the center')

		string=self.modelPrint()
		fig.text(0.0,0.14,string+'Av='+str(self.extinction)+'\n'+'dist='+str(self.dist/pc)+'\n',color='r')
		fig.savefig(modelname+'.png', bbox_inches='tight',dpi=300)

		if show:
			plt.show()
			
	def plotSim(self,dist_pc=140,inc=3,extinction=0,show=False):
		self.dist=dist_pc*pc
		self.inc=inc
		self.extinction=extinction
		modelname = self.folder+self.name
		self.mo = ModelOutput(modelname+'.rtout')

		#tracy_dust = np.loadtxt('Tracy_models/OH5.par')
		#chi = np.loadtxt('kmh94_3.1_full.chi')
		#wav = np.loadtxt('kmh94_3.1_full.wav')
		#Chi = interp1d(wav,chi,kind='linear')



		fig = plt.figure(figsize=(20,14))
		ax=fig.add_subplot(1,3,1)
		sed = self.mo.get_sed(aperture=-1, inclination='all', distance=self.dist)
		#print tracy_dust[11,1],Cext(sed.wav[-1]),Cext(sed.wav[-1])/tracy_dust[11,1]
		#tau = self.extinction*Chi(sed.wav)/Chi(0.550)/1.086
		#print Cext(sed.wav)/tracy_dust[11,1]
		#ext = np.array([np.exp(-tau) for i in range(sed.val.shape[0])])
		#print tau,np.exp(-tau)
		ax.loglog(sed.wav, sed.val.transpose(), color='black')
		ax.set_title(modelname+'_seds, Av='+str(self.extinction))
		ax.set_xlim(sed.wav.min(), 1300)
		ax.set_ylim(1e-13, 1e-7)
		ax.set_xlabel(r'$\lambda$ [$\mu$m]')
		ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
		#self.plotData(ax,sourcename)
		ax.set_xscale('log')
		ax.set_yscale('log')

		#ax.set_ylabel(r'$F_{Jy}$ [Jy]')
		#plt.legend(loc=4)

#		ax=fig.add_subplot(2,3,2)
#		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist)
#		ext=np.exp(-tau)
#		ax.loglog(sed.wav, sed.val.transpose()*ext.T, lw=3,color='black',label='source_total')
#		ax.set_xlim(sed.wav.min(), 1300)
#		ax.set_ylim(1e-13, 1e-7)  ### for lamFlam
#		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='source_emit')
#		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='blue',label='source_emit')
#		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='source_scat')
#		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='teal',label='source_scat')
#		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='dust_emit')
#		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='red',label='dust_emit')
#		sed = self.mo.get_sed(aperture=-1, inclination=self.inc, distance=self.dist,component='dust_scat')
#		ax.loglog(sed.wav, sed.val.transpose()*ext.T, color='orange',label='dust_scat')
#		#self.plotData(ax,sourcename)
#		ax.set_xscale('log')
#		ax.set_yscale('log')
#		ax.set_title('seds_inc=inc')
#		ax.set_xlabel(r'$\lambda$ [$\mu$m]')
#		ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
#		#ax.set_ylabel(r'$F_{Jy}$ [Jy]')
#		leg = ax.legend(loc=4,fontsize='small')
		#leg = plt.gca().get_legend()
		#plt.setp(leg.get_text(),fontsize='small')
		# Extract the quantities
		g = self.mo.get_quantities()
		
		# Get the wall positions for r and theta
		rw, tw = g.r_wall / au, g.t_wall

		# Make a 2-d grid of the wall positions (used by pcolormesh)
		R, T = np.meshgrid(rw, tw)

		# Calculate the position of the cell walls in cartesian coordinates
		X, Z = R * np.sin(T), R * np.cos(T)

		# Make a plot in (x, z) space for different zooms
		from matplotlib.colors import LogNorm,PowerNorm
		# Make a plot in (r, theta) space
		ax = fig.add_subplot(1, 3, 2)
		if g.shape[-1]==2:
			c = ax.pcolormesh(X, Z, g['temperature'][0].array[0, :, :]+g['temperature'][1].array[0, :, :],norm=PowerNorm(gamma=0.5,vmin=1,vmax=500))
		else :
			c = ax.pcolormesh(X, Z, g['temperature'][0].array[0, :, :],norm=PowerNorm(gamma=0.5,vmin=1,vmax=500))
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlim(X.min(), X.max())
		ax.set_ylim(Z.min(), Z.max())
		ax.set_xlabel('x (au)')
		ax.set_ylabel('z (au)')
		#ax.set_yticks([np.pi, np.pi * 0.75, np.pi * 0.5, np.pi * 0.25, 0.])
		#ax.set_yticklabels([r'$\pi$', r'$3\pi/4$', r'$\pi/2$', r'$\pi/4$', r'$0$'])
		cb = fig.colorbar(c)
		ax.set_title('Temperature structure')
		cb.set_label('Temperature (K)')
		#fig.savefig(modelname+'_temperature_spherical_rt.png', bbox_inches='tight')


		ax = fig.add_subplot(1, 3, 3)
		if g.shape[-1]==2:
			c = ax.pcolormesh(X, Z, g['density'][0].array[0, :, :]+g['density'][1].array[0, :, :],norm=LogNorm(vmin=1e-22,vmax=g['density'][0].array[0, :, :].max()))
		else :
			c = ax.pcolormesh(X, Z, g['density'][0].array[0, :, :],norm=LogNorm(vmin=1e-22,vmax=g['density'][0].array[0, :, :].max()))
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlim(X.min(), X.max())
		ax.set_ylim(Z.min(), Z.max())
		ax.set_xlabel('x (au)')
		ax.set_ylabel('z (au)')
		ax.set_title('Density structure')
		cb = fig.colorbar(c)
		cb.set_label('Density (g/cm2)')

#		### plot the convolved image with the 37 micron filter (manually set to slice 18 of the cube - this would change with wavelength coverage)
#		ax = fig.add_subplot(2, 3, 5)
#		self.image = self.mo.get_image(inclination=inc,distance=self.dist,units='Jy')
#		fits.writeto(modelname+'_inc_'+str(inc)+'.fits',self.image.val.swapaxes(0,2).swapaxes(1,2),clobber=True)

#		### need to convolve the image with a Gaussian PSF
#		pix = 2.*self.limval/au/self.Npix # in AU/pix
#		pix_asec = pix/(self.dist/pc) # in asec/pix
#		airy_asec = 3.5 #asec
#		airy_pix = airy_asec/pix_asec # in pix
#		gauss_pix = airy_pix/2.35 # in Gaussian std 
#		print "Gaussian std: ",gauss_pix

#		from scipy.ndimage.filters import gaussian_filter as gauss
#		#print [(i,sed.wav[i]) for i in range(len(sed.wav))]

#		img37 = self.image.val[:,:,18]
#		convol = gauss(img37,gauss_pix,mode='constant',cval=0.0)
#		Nc = self.Npix/2
#		hw = min(int(20./pix_asec),Nc) #(max is Nc)
#		#ax.imshow(img37,norm=LogNorm(vmin=1e-20,vmax=img37.max()))
#		#ax.imshow(img37,interpolation='nearest')
#		#ax.imshow(convol,norm=LogNorm(vmin=1e-20,vmax=img37.max()))
#		#ax.imshow(convol,interpolation='nearest',norm=LogNorm(vmin=1e-20,vmax=img37.max()))
#		ax.imshow(convol[Nc-hw:Nc+hw,Nc-hw:Nc+hw],interpolation='nearest',origin='lower',cmap=plt.get_cmap('gray'))
#		airy_disk = plt.Circle((airy_pix*1.3,airy_pix*1.3),airy_pix,color=colors[3])		
#		ax.add_artist(airy_disk)
#		ax.text(airy_pix*3,airy_pix*1.3/2.0,'SOFIA 37um Airy disk',color=colors[3])
#		ax.set_title('Convolved image')
#		fits.writeto(modelname+'_inc_'+str(inc)+'_convol37.fits',convol,clobber=True)

#		### draw a cross-section of the image to show the spatial extension in linear scale, to compare with what we observe in the model.
#		ax = fig.add_subplot(2, 3, 6)
#		ax.plot(range(Nc-hw,Nc+hw),convol[Nc-hw:Nc+hw,Nc-1],label='cross-section 1')
#		ax.plot(range(Nc-hw,Nc+hw),convol[Nc-1,Nc-hw:Nc+hw],label='cross-section 2')
#		maxconvol = convol[Nc-hw:Nc+hw,Nc-1].max()
#		gauss = np.exp( -(np.array(range(-hw,hw))**2 / (2. * gauss_pix**2)))
#		gauss/= gauss.max()
#		gauss*=maxconvol
#		ax.plot(range(Nc-hw,Nc+hw),gauss,label='SOFIA beam')
#		leg = ax.legend(loc=2,fontsize='small')
#		#leg = plt.gca().get_legend()
#		#plt.setp(leg.get_text(),fontsize='small')
#		ax.set_title('Cross section at the center')

		string=self.modelPrint()
		fig.text(0.0,0.14,string+'Av='+str(self.extinction)+'\n'+'dist='+str(self.dist/pc)+'\n',color='r')
		fig.savefig(modelname+'.png', bbox_inches='tight',dpi=300)

		if show:
			plt.show()

def modelLoad(folder,name):
	YSOModel = pickle.load(open(folder+name+'.mod','r'))
	return YSOModel

def modelDump(YSOmodel):
	sp.call('rm %s.mod ' % (YSOmodel.folder+YSOmodel.name),shell=True)
	time.sleep(2)
	pickle.dump(YSOmodel,open(YSOmodel.folder+YSOmodel.name+'.mod','wb'))
	time.sleep(2)
