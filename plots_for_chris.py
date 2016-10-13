#Will need to compile things together.
import os
import astropy.cosmology
import multiprocessing
import astropy as ast
import pylab as pl
from astropy.io import fits
from matplotlib.colors import LogNorm
from numpy.ma import MaskedArray
#import aplpy as ap
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import csv
import re
import urllib
from astropy.coordinates import SkyCoord
import matplotlib.image as Image
from astropy import units as u
import os,sys
from tempfile import mkstemp
from shutil import move
from os import remove, close
import time
from scipy import stats
from scipy.stats import ks_2samp
from scipy.integrate import quad
from astropy.modeling.models import Sersic1D
from scipy.stats import kstest
import astropy.stats as asts
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import matplotlib.mlab as mlab
import math
from numpy import arccos
from numpy import cos
from numpy import sin
from astropy.cosmology import z_at_value
from matplotlib.pyplot import cm
import cosmolopy as cp
import sklearn as sk
from scipy.stats import spearmanr
################################################################################
#Functions
################################################################################
def csvsave(filenameandurl,header,file_input):
	x = 0
	myfile = open(np.str(filenameandurl), 'wb')
	wr = csv.writer(myfile, quoting=csv.QUOTE_NONE, delimiter=',')
	for g in range(0,len(file_input)+1):
		if g ==0:
			wr.writerow(header)
		elif x < len(file_input):
			wr.writerow(file_input[x])
			x = x + 1
	return ;

def intstr(input):
	output = np.str(np.int(input))
	return output;

def getheadparam(file_id,extension,param):
	data = fits.open(file_id)
	variable = data[extension].header[param]
	return variable;
# def get_sdss_image(urls,band):

def get_sdss_jpeg(rerun,run,camcol,field,location):
	rerun = intstr(rerun)
	camcol = intstr(camcol)
	field = intstr(field)
	run = intstr(run)
	if np.int(field) < 10:
		getfile = ['http://data.sdss3.org/sas/dr12/boss/photoObj/frames/'+rerun+'/'+run+'/'+camcol+'/frame-irg-00'+run+'-'+camcol+'-000'+field+'.jpg']
	elif np.int(field) >10 and np.int(field) <100:
		getfile = ['http://data.sdss3.org/sas/dr12/boss/photoObj/frames/'+rerun+'/'+run+'/'+camcol+'/frame-irg-00'+run+'-'+camcol+'-00'+field+'.jpg']
	else:
		getfile = ['http://data.sdss3.org/sas/dr12/boss/photoObj/frames/'+rerun+'/'+run+'/'+camcol+'/frame-irg-00'+run+'-'+camcol+'-0'+field+'.jpg']
	urllib.urlretrieve(getfile[0], filename=location)
	return ;

def join_table(ra1,dec1,z1,ra2,dec2,id2): #note at least one needs object id
	catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
	match = []
	for i in range(0,len(ra1)):
		c = SkyCoord(ra=ra1[i]*u.degree, dec=dec1[i]*u.degree)
		idx1, d2d1, d3d1 = c.match_to_catalog_sky(catalog)
		if d2d1.degree < 0.5: #thirty arcsecond match
			match.append([ra1[i],dec1[i],z1[i],ra2[idx1],dec2[idx1],id2[idx1],np.float(d2d1.degree)])
			print intstr(i) # this will provide the matching index for the catalogue. Then , you can use these to pick out the clusters and group similar LRGs
	return match ; #goes for the nearest match

def extract_fits_table(file_location):
	datafile = fits.open(file_location)
	datafile = datafile[1].data
	return datafile ;

def arc_to_kpc(z,r): #z and r must be the same length
	if len(z) == 1:
		arc_to_kpc = cosmo.kpc_proper_per_arcmin(z)
		r_kpc = arc_to_kpc*(r/60)
	elif len(z) > 0 and len(z) != 1:
		for i in range(0,len(z)):
			arc_to_kpc = cosmo.kpc_proper_per_arcmin(z[i])
			r_kpc = arc_to_kpc*(r[i]/60)
	else:
		r_kpc == np.nan
	return r_kpc ;

def launch_DS9_multiframe(image): #MUST BE IN IMAGE DIRECTORY
	os.system("ds9 -multiframe "+image)
	return ;

def mask_crunch2(array_1, array_2, mask1, mask2):
	array1 = np.array(array_1[mask1]).tolist()
	array2 = np.array(array_2[mask2]).tolist()
	total = array1 + array2
	return total;

def mask_crunch(array, mask1, mask2):
	array1 = np.array(array[mask1]).tolist()
	array2 = np.array(array[mask2]).tolist()
	total = array1 + array2
	return total;
################################################################################
#v disp - comp with lit for guo and zhao, z as well (z_diff binned, guo and zhao) - err?

#example of parameter space? explain what I'm doing? Plot Re vs n with errors, three mag bins
#explain that she's running sims at the minute to explore capabilities of software

#See if can find another way of separating morphs - zhao looked at RFF

#Redshift range of sample - working on calibrating mags for newest set

#mass = 1.15 + 0.70*(kcor_mg -  kcor_mi) - 0.4*(k_cor_Mi) - for errors in mass
################################################################################
data = extract_fits_table('bcg_final_v2_dust_zbwt_only.fits')

scale_type = data.field('SCALE_TYPE')
stellar_masses = data.field('log_BCG_Mass_Taylor_dust')#previously, the 'v' stood for de Vaucouleurs
m200_clus_gapper = data.field('M200_SIGMA_GAPPER')
m200_gapper_err = data.field('GAPPER_SIGMA')
m200_clus_bwt = data.field('M200_SIGMA_BIWEIGHT')
m200_bwt_err = data.field('BIWEIGHT_SIGMA')
z_clus = data.field('ZBWT')
Mi_kcorr_dust = data.field('Mi_s_dust_2a')
mg_kcorr_dust = data.field('mg_s_dust_2')
mi_kcorr_dust = data.field('mi_s_dust_2')

sigma_flag = data.field('SIG_FLAG')
sigma_bwt = [sigma_flag == 'BIWEIGHT']
sigma_gapper = [sigma_flag == 'GAPPER']

stellar_mass_total = np.array(mask_crunch(stellar_masses,sigma_gapper,sigma_bwt))
z_total = mask_crunch(z_clus,sigma_gapper,sigma_bwt)
m200_total = np.array(mask_crunch2(m200_clus_gapper,m200_clus_bwt,sigma_gapper,sigma_bwt))
m200_err_total = mask_crunch2(m200_gapper_err,m200_bwt_err,sigma_gapper,sigma_bwt)
m200_err_total = np.array((1.0/0.7)*(1.2*10**15)*((np.array(m200_err_total)/1000.0)**3)*((np.sqrt(0.7 + 0.3*(1 + np.array(z_total))**3.0))**(-1.0)))

mg_kcorr_dust_total = mask_crunch(mg_kcorr_dust,sigma_gapper,sigma_bwt)
mi_kcorr_dust_total = mask_crunch(mi_kcorr_dust,sigma_gapper,sigma_bwt)
Mi_kcorr_dust_total = mask_crunch(Mi_kcorr_dust,sigma_gapper,sigma_bwt)
#stellar_mass_errors = 1.15 + 0.70*(mg_kcorr_dust_total -  mi_kcorr_dust_total) - 0.4*(Mi_kcorr_dust_total) #- for errors in mass
m200_model = np.array([12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0])
m, b = np.polyfit(np.log10(m200_total), stellar_mass_total, 1)
ax = plt.plot(np.log10(m200_total),stellar_mass_total,'k.',markersize=6)
plt.plot(np.log10(m200_total),stellar_mass_total,'k.',markersize=6)
plt.errorbar(np.log10(m200_total),stellar_mass_total,xerr=0.434*(np.array(m200_err_total)/np.array(m200_total)),fmt='k.')
plt.plot(m200_model, m*m200_model + b, 'b--',label=r'($\alpha$,  $\beta $) = '+'('+np.str(np.round(m,2))+', '+np.str(np.round(b,2))+')',linewidth=2)
s, p = spearmanr(np.log10(m200_total),stellar_mass_total)
plt.annotate(r'$\mathrm{r_{s}}$ = '+np.str(np.round(s,3)), xy=(0.2, 0.84), xycoords='axes fraction', fontsize=14, horizontalalignment='right', verticalalignment='bottom')
plt.annotate('log10(p) = '+np.str(np.round(np.log10(p),3)), xy=(0.32, 0.78), xycoords='axes fraction', fontsize=14, horizontalalignment='right', verticalalignment='bottom')
plt.legend(loc='upper left',frameon=False)
plt.xlabel('log(Cluster Mass)',fontsize=14)
plt.ylabel('log(BCG Mass)',fontsize=14)
plt.show()
