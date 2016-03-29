'''
This script uses a simple Chi Squared estimate to estimate which model from the grid best matches the data from an SED
'''
import photometry as p
import pickle
import numpy as np
from astropy.table import Table
import os
from scipy.interpolate import interp1d,interp2d,SmoothBivariateSpline
from hyperion.model import ModelOutput
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import seaborn as sns

sns.set(context='talk',style='ticks')
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.labelsize']=16
mpl.rcParams['legend.fontsize']=14
mpl.rcParams['font.size']=14


# set the target
target_list = ['IRAS20050.1','IRAS20050.2','IRAS20050.3','IRAS20050.4','IRAS20050.5','IRAS20050.6','IRAS20050.7','NGC2071.1','NGC2071.2','NGC2071.3','NGC2071.4','NGC2071.5']
#'IRAS20050.1','IRAS20050.2','IRAS20050.3','IRAS20050.4','IRAS20050.5','IRAS20050.6','IRAS20050.7',
def distance(target):
	if "IRAS20050" in target:
		return 700.*pc
	elif "NGC2071" in target:
		return 390.*pc
		
red = sns.xkcd_rgb['pale red']
green = sns.xkcd_rgb['medium green']
blue = sns.xkcd_rgb['denim blue']
amber = sns.xkcd_rgb['amber']

# load the data points for that target
folder_export="/n/a2/mrizzo/Dropbox/SOFIA/Processed_Data/"
sourcetable = pickle.load(open(folder_export+"totsourcetable_fits.data","r"))

# Process data points
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
#HERSCHEL = plots(['H70','H160','H250','H350','H500'],[70,160,250,350,500],colors[7],markers[7],'HERSCHEL')


# set up extinction
extinctions = range(30)
d = SphericalDust()
d.read('d03_5.5_3.0_A.hdf5')
chi = d.optical_properties.chi
chi = chi[::-1]
wav = d.optical_properties.wav
wav = wav[::-1]
Chi = interp1d(wav,chi,kind='linear')

#inclinations = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]
angles=np.arccos(np.linspace(0,1.,20))*180./np.pi
inclinations=angles[::-1]


# set the wavelengths
names = TwoMASS+Spitzer+SOFIA
wl=wlTwoMASS+wlSpitzer+wlSOFIA
wl_table = Table(names = names,dtype=['f8' for col in wl])
wl_table.add_row(wl)
N = len(wl)


# load the grid
name = ['model']
folder = ['Gridfinal/']
filename = folder[0]+name[0]+"_eval.grid.dat"
if os.path.exists(filename):
	f = open(filename,'r')
	grid_tot = pickle.load(f)
	f.close()

else:
	#TODO: do something when the file isn't found!
	print "No grid reference file found!"
	
grid_ext = grid_tot.group_by('ext')
ext = 0.0


# select only the models taken at zero extinction
for key,gridi in zip(grid_ext.groups.keys,grid_ext.groups):
	if key['ext']==ext:
		#grid = gridi
		gridi_mdisk = gridi.group_by('disk_mass')
		for keym,gridm in zip(gridi_mdisk.groups.keys,gridi_mdisk.groups):
			if keym['disk_mass']==0.001:
				grid=gridm
	
# compute color for all models, extinctions and inclinations
columnlist = ['j','h','ks','w1','i1','i2','w2','i3','i4','F11','w3','Fnu_12','F19','w4','m1','Fnu_25','F31','F37']
grid['i4-i2'] = grid['i4'] - grid['i2']
grid['F37-F19'] = grid['F37'] - grid['F19']
grid['F31-i3'] = grid['F31'] - grid['i3']
sourcetable = sourcetable[['SOFIA_name']+columnlist]
sources = sourcetable.group_by('SOFIA_name')

columns = ['name', 'folder','T','M_sun','env_rmax','env_rmin','disk','disk_mass','disk_rmax',\
					'cav_theta','beta','L_sun','env_mass','env_power']
grid.sort('i4-i2')
#grid[columns+['inc','ext','i4-i2','F37-F19']].more()
# Main plotting function
def plot_color_color(variable,ext,grid=grid,label=""):
	fig = plt.figure(figsize=(25,15))
	ax = fig.add_axes([0.08,0.1,0.9,0.8])
	ax.grid(True)
	color_code = np.array(grid[variable])
	colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
	ax.set_axis_bgcolor((0.8,0.8,0.8))
	xmin = min(grid['i4-i2'])
	xmax = max(grid['i4-i2'])
	ymin = min(grid['F37-F19'])
	ymax = max(grid['F37-F19'])
	ax.set_xlabel('[i4]-[i2]')
	ax.set_ylabel('[F37]-[F19]')
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.scatter(grid['i4-i2'],grid['F37-F19'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')
	ax2 = fig.add_axes([0.1,0.15,0.04,0.4])
	norm = mpl.colors.Normalize(vmin=min(color_code),vmax = max(color_code))
	cbar = mpl.colorbar.ColorbarBase(ax2,cmap = "Reds_r",norm=norm)
	cbar.set_label(variable)

	for target in target_list:
		print target
		for key,source in zip(sources.groups.keys,sources.groups):
			if target == source['SOFIA_name'][0]:	
				print "found!"
				datapoints = source[names]
				color1 = np.log10(datapoints['i4'][0]) - np.log10(datapoints['i2'][0]) #use the modified expression here?
				color2 = np.log10(datapoints['F37'][0]) - np.log10(datapoints['F19'][0])
				if "NGC2071" in target:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=blue,edgecolor='k',marker='^')
				else:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=green,edgecolor='k',marker='^')
	fig.suptitle("Extinction="+str(ext)+", "+variable+"_"+label)
	fig.savefig(folder[0]+variable+"_ext="+str(ext)+label+".png",dpi=300)

# Mosaic plotting function
def plot_color_color_mosaic(ax,variable,ext,grid=grid,label=""):
	#fig = plt.figure(figsize=(25,15))
	#ax = fig.add_axes([0.08,0.1,0.9,0.8])
	ax.grid(True)
	if variable=='inc':
		color_code = np.cos(np.array(grid[variable])*np.pi/180.)
	else:
		color_code = np.array(grid[variable])
	colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
	ax.set_axis_bgcolor((0.8,0.8,0.8))
	xmin = min(grid['i4-i2'])
	xmax = max(grid['i4-i2'])
	ymin = min(grid['F37-F19'])
	ymax = max(grid['F37-F19'])
	ax.set_xlabel('[i4]-[i2]')
	ax.set_ylabel('[F37]-[F19]')
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.scatter(grid['i4-i2'],grid['F37-F19'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')
	ax.set_title(variable+"_ext="+str(ext))
	pos = ax.get_position()
	ax2 = fig.add_axes([pos.x0+0.1*pos.width,pos.y0+0.15*pos.height,0.04*pos.width,0.4*pos.height])
	norm = mpl.colors.Normalize(vmin=min(color_code),vmax = max(color_code))
	cbar = mpl.colorbar.ColorbarBase(ax2,cmap = "Reds_r",norm=norm)
	cbar.set_label(variable)

	for target in target_list:
		print target
		for key,source in zip(sources.groups.keys,sources.groups):
			if target == source['SOFIA_name'][0]:	
				print "found!"
				datapoints = source[names]
				color1 = np.log10(datapoints['i4'][0]) - np.log10(datapoints['i2'][0]) #use the modified expression here?
				color2 = np.log10(datapoints['F37'][0]) - np.log10(datapoints['F19'][0])
				if "NGC2071" in target:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=blue,edgecolor='k',marker='^')
				else:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=green,edgecolor='k',marker='^')
	#fig.suptitle("Extinction="+str(ext)+", "+variable)
	#fig.savefig(folder[0]+variable+"_ext="+str(ext)+".png",dpi=300)

def plot_color_color_mosaic3D(ax,variable,ext,grid=grid,label=""):
	#fig = plt.figure(figsize=(25,15))
	#ax = fig.add_axes([0.08,0.1,0.9,0.8])
	ax.grid(True)
	if variable=='inc':
		color_code = np.cos(np.array(grid[variable])*np.pi/180.)
	else:
		color_code = np.array(grid[variable])
	colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
	ax.set_axis_bgcolor((0.8,0.8,0.8))
	xmin = min(grid['i4-i2'])
	xmax = max(grid['i4-i2'])
	ymin = min(grid['F37-F19'])
	ymax = max(grid['F37-F19'])
	zmin = min(grid['F31-i3'])
	zmax = max(grid['F31-i3'])
	ax.set_xlabel('[i4]-[i2]')
	ax.set_ylabel('[F37]-[F19]')
	ax.set_zlabel('[F31]-[i3]')
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.set_zlim([zmin,zmax])
	ax.scatter(grid['i4-i2'],grid['F37-F19'],grid['F31-i3'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')
	ax.set_title(variable+"_ext="+str(ext)+label)
	pos = ax.get_position()
	ax2 = fig.add_axes([pos.x0+0.1*pos.width,pos.y0+0.15*pos.height,0.04*pos.width,0.4*pos.height])
	norm = mpl.colors.Normalize(vmin=min(color_code),vmax = max(color_code))
	cbar = mpl.colorbar.ColorbarBase(ax2,cmap = "Reds_r",norm=norm)
	cbar.set_label(variable)

	for target in target_list:
		print target
		for key,source in zip(sources.groups.keys,sources.groups):
			if target == source['SOFIA_name'][0]:	
				print "found!"
				datapoints = source[names]
				color1 = np.log10(datapoints['i4'][0]) - np.log10(datapoints['i2'][0]) #use the modified expression here?
				color2 = np.log10(datapoints['F37'][0]) - np.log10(datapoints['F19'][0])
				color3 = np.log10(datapoints['F31'][0]) - np.log10(datapoints['i3'][0])
				if "NGC2071" in target:
					ax.scatter(color1,color2,color3,s=200,alpha=0.9,color=blue,edgecolor='k',marker='^')
				else:
					ax.scatter(color1,color2,color3,s=200,alpha=0.9,color=green,edgecolor='k',marker='^')

from mpl_toolkits.mplot3d import Axes3D

#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')
#variable = 'env_mass'
#color_code = np.array(grid[variable])
#colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
#ax.scatter(grid['i4-i2'],grid['F37-F19'],grid['F31-i3'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')

### individual plots
## variable = envelope mass
#plot_color_color('env_mass',ext)
## variable = luminosity
#plot_color_color('L_sun',ext)
## variable = envelope power
#plot_color_color('env_power',ext)
## variable = inclination
#plot_color_color('inc',ext)

### 2D mosaics
#fig = plt.figure(figsize=(25,15))
#gs = gridspec.GridSpec(2,2)
#gs.update(left=0.05,right=0.95,top=0.97,bottom=0.06)
## variable = envelope mass
#plot_color_color_mosaic(plt.subplot(gs[0]),'env_mass',ext)
## variable = luminosity
#plot_color_color_mosaic(plt.subplot(gs[1]),'L_sun',ext)
## variable = envelope power
#plot_color_color_mosaic(plt.subplot(gs[2]),'env_power',ext)
## variable = inclination
#plot_color_color_mosaic(plt.subplot(gs[3]),'inc',ext)

### 3D mosaics
#fig = plt.figure(figsize=(25,15))
#gs = gridspec.GridSpec(2,2)
#gs.update(left=0.05,right=0.95,top=0.97,bottom=0.06)
## variable = envelope mass
#plot_color_color_mosaic3D(plt.subplot(gs[0],projection='3d'),'env_mass',ext)
## variable = luminosity
#plot_color_color_mosaic3D(plt.subplot(gs[1],projection='3d'),'L_sun',ext)
## variable = envelope power
#plot_color_color_mosaic3D(plt.subplot(gs[2],projection='3d'),'env_power',ext)
## variable = inclination
#plot_color_color_mosaic3D(plt.subplot(gs[3],projection='3d'),'inc',ext)


#fig.savefig(folder[0]+"mosaic_ext="+str(ext)+".png",dpi=300)


# Plotting function for other color-color diagram
def plot_color_color2(variable,ext,grid=grid,label="c2"):
	fig = plt.figure(figsize=(25,15))
	ax = fig.add_axes([0.08,0.1,0.9,0.8])
	ax.grid(True)
	color_code = np.array(grid[variable])
	colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
	ax.set_axis_bgcolor((0.8,0.8,0.8))
	xmin = min(grid['F37-F19'])
	xmax = max(grid['F37-F19'])
	ymin = min(grid['F31-i3'])
	ymax = max(grid['F31-i3'])
	ax.set_xlabel('[F37]-[F19]')
	ax.set_ylabel('[F31]-[i3]')
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.scatter(grid['F37-F19'],grid['F31-i3'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')
	ax2 = fig.add_axes([0.5,0.15,0.3,0.04])
	norm = mpl.colors.Normalize(vmin=min(color_code),vmax = max(color_code))
	cbar = mpl.colorbar.ColorbarBase(ax2,cmap = "Reds_r",norm=norm,orientation='horizontal')
	cbar.set_label(variable)

	for target in target_list:
		print target
		for key,source in zip(sources.groups.keys,sources.groups):
			if target == source['SOFIA_name'][0]:	
				print "found!"
				datapoints = source[names]
				color1 = np.log10(datapoints['F37'][0]) - np.log10(datapoints['F19'][0]) #use the modified expression here?
				color2 = np.log10(datapoints['F31'][0]) - np.log10(datapoints['i3'][0])
				if "NGC2071" in target:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=blue,edgecolor='k',marker='^')
				else:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=green,edgecolor='k',marker='^')
	fig.suptitle("Extinction="+str(ext)+", "+variable+"_"+label)
	fig.savefig(folder[0]+variable+"_ext="+str(ext)+label+".png",dpi=300)
	
# Plotting function for other color-color diagram
def plot_color_color3(variable,ext,grid=grid,label="c3"):
	fig = plt.figure(figsize=(25,15))
	ax = fig.add_axes([0.08,0.1,0.9,0.8])
	ax.grid(True)
	color_code = np.array(grid[variable])
	colors = (color_code-min(color_code))/(max(color_code)-min(color_code))
	ax.set_axis_bgcolor((0.8,0.8,0.8))
	xmin = min(grid['i4-i2'])
	xmax = max(grid['i4-i2'])
	ymin = min(grid['F31-i3'])
	ymax = max(grid['F31-i3'])
	ax.set_xlabel('[i4]-[i2]')
	ax.set_ylabel('[F31]-[i3]')
	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.scatter(grid['i4-i2'],grid['F31-i3'],cmap="Reds_r",s=20,alpha=0.5,c=colors,edgecolor='none')
	ax2 = fig.add_axes([0.5,0.15,0.3,0.04])
	norm = mpl.colors.Normalize(vmin=min(color_code),vmax = max(color_code))
	cbar = mpl.colorbar.ColorbarBase(ax2,cmap = "Reds_r",norm=norm,orientation='horizontal')
	cbar.set_label(variable)

	for target in target_list:
		print target
		for key,source in zip(sources.groups.keys,sources.groups):
			if target == source['SOFIA_name'][0]:	
				print "found!"
				datapoints = source[names]
				color1 = np.log10(datapoints['i4'][0]) - np.log10(datapoints['i2'][0]) #use the modified expression here?
				color2 = np.log10(datapoints['F31'][0]) - np.log10(datapoints['i3'][0])
				if "NGC2071" in target:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=blue,edgecolor='k',marker='^')
				else:
					ax.scatter(color1,color2,s=200,alpha=0.9,color=green,edgecolor='k',marker='^')
	fig.suptitle("Extinction="+str(ext)+", "+variable+"_"+label)
	fig.savefig(folder[0]+variable+"_ext="+str(ext)+label+".png",dpi=300)


## variable = envelope mass
#var = 'env_mass'
#plot_color_color(var,ext)
#plot_color_color2(var,ext)
#plot_color_color3(var,ext)
## variable = envelope mass
#var = 'L_sun'
#plot_color_color(var,ext)
#plot_color_color2(var,ext)
#plot_color_color3(var,ext)

grid.write(folder[0]+name[0]+"_test_export.grid.txt",format='ascii.fixed_width')

# select only the models taken at zero extinction
grid_inc = grid.group_by('inc')
for keyi,gridi in zip(grid_inc.groups.keys,grid_inc.groups):
	print keyi['inc']
	print grid_inc.groups.keys['inc']
	if keyi['inc'] == grid_inc.groups.keys['inc'][0] or keyi['inc'] == grid_inc.groups.keys['inc'][-1] or keyi['inc'] == grid_inc.groups.keys['inc'][len(grid_inc.groups.keys)/2]:
		# variable = envelope mass
		print "current:",keyi['inc']
		var = 'L_sun'
		plot_color_color(var,ext,grid=gridi,label="inc_"+str(keyi['inc']))
		plot_color_color2(var,ext,grid=gridi,label="inc_"+str(keyi['inc'])+"c2")
		plot_color_color3(var,ext,grid=gridi,label="inc_"+str(keyi['inc'])+"c3")
		fig = plt.figure(figsize=(25,15))
		gs = gridspec.GridSpec(1,1)
		gs.update(left=0.05,right=0.95,top=0.97,bottom=0.06)
		# variable = envelope mass
		plot_color_color_mosaic3D(plt.subplot(gs[0],projection='3d'),var,ext,grid=gridi,label="inc_"+str(keyi['inc']))
		fig.savefig(folder[0]+"mosaic_ext="+str(ext)+"_inc_"+str(keyi['inc'])+".png",dpi=300)

		


plt.show()



