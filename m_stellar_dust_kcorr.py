import os
import astropy.cosmology
import multiprocessing
import astropy as ast
import pylab as pl
from astropy.io import fits
from matplotlib.colors import LogNorm
from numpy.ma import MaskedArray
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import EllipticalAperture
#import aplpy as ap
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import csv
import re
import urllib
from astropy.coordinates import SkyCoord
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
import photutils
from astropy.modeling import models, fitting
from scipy.stats import spearmanr
from matplotlib import patches
import matplotlib.pyplot as plt
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from calc_kcor import calc_kcor #make sure this is in your path before running it
################################################################################
def doFitting(FileList):
	command = ("./galfit -o3 "+FileList)
	print command
	os.system(command)

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
def getheadparam(file_id,extension,param):
	data = fits.open(file_id)
	variable = data[extension].header[param]
	return variable;

def intstr(input):
	output = np.str(np.int(input))
	return output;

def replace(file_path, new_path, pattern, subst):
	#Create temp file
	fh, abs_path = mkstemp()
	with open(abs_path,'w') as new_file:
		with open(file_path) as old_file:
			for line in old_file:
				new_file.write(line.replace(pattern, subst))
	close(fh)
	move(abs_path, new_path)

def kcorr_dust_mass(g,r,i,z_bcg,ebv):
	#SaF 2011 coeffcients
	g_band_ex = 3.303
	r_band_ex = 2.285
	i_band_ex = 1.698
	#dust correction
	Ag = g_band_ex*ebv
	Ar = r_band_ex*ebv
	Ai = i_band_ex*ebv
	#Dust corrected magnitudes
	gmag_dust = g - Ag
	rmag_dust = r - Ar
	imag_dust = i - Ai
	#colours
	g_i = gmag_dust - imag_dust
	g_r = gmag_dust - rmag_dust
	#k-corrected magnitudes
	lum_mpc = cosmo.luminosity_distance(z_bcg) #bcg z
	lum_mpc = lum_mpc.value
	Mg = gmag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
	Mr = rmag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
	Mi = imag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
	kcor_r = calc_kcor('r',z_bcg,'g - r',g_r)
	kcor_g = calc_kcor('g',z_bcg,'g - r',g_r)
	kcor_i = calc_kcor('i',z_bcg,'g - i',g_i)
	k_cor_Mr = Mr - kcor_r
	k_cor_Mi = Mi - kcor_i
	k_cor_Mg = Mg - kcor_g
	kcor_mr = rmag_dust - kcor_r
	kcor_mg = gmag_dust - kcor_g
	kcor_mi = imag_dust - kcor_i
	mass = 1.15 + 0.70*(kcor_mg -  kcor_mi) - 0.4*(k_cor_Mi) #restframe colours. Using Taylor et al., 2015 relation rather than 2011
	return Ag,Ar,Ai,kcor_g,kcor_r,kcor_i,kcor_mg,kcor_mr,kcor_mi,mass ;

################################################################################
#K-corrections - first shot -  calc_kcor(filter_name, redshift, colour_name, colour_value):
#reddening calculator - do this BEFORE you k-correct
url = '/home/kate/Desktop/'
dust_values = fits.open(url+'complete_run_full_sample_same_scale_nb.fits')
dust_values = dust_values[1].data
g_band_ex = 3.303 #Schlafly et al., 2011 bandpass conversions (Rv ~ 3.1)
r_band_ex = 2.285
i_band_ex = 1.698

refname = dust_values.field('CLUS_ID_2')
E_B_V = dust_values.field('E_B_V_SandF')
g_mags = dust_values.field('GAL_sdss_g_modS_C1_MAG')
r_mags = dust_values.field('GAL_sdss_r_modS_C1_MAG')
i_mags = dust_values.field('GAL_sdss_i_modS_C1_MAG')
z_clus = dust_values.field('Z_BWT_1')

#PARAMAETERS
Ag_out = []
Ar_out = []
Ai_out = []
k_cor_mr_out = []
k_corr_mg_out = []
k_corr_mi_out = []
k_corr_mr_out = []
k_corr_Mi_out = []
k_corr_Mg_out = []
k_corr_Mr_out = []
mass_out = []
kcor_r_out = []
kcor_i_out = []
kcor_g_out = []
ref_id = []
for i in range(0,len(refname)):
	if g_mags[i] != -9999 and r_mags[i] != -9999 and i_mags[i] != -9999:
		#REDDENING Coefficients,
		Ag = g_band_ex*E_B_V[i]
		Ar = r_band_ex*E_B_V[i]
		Ai = i_band_ex*E_B_V[i]
		#Dust corrected magnitudes
		gmag_dust = g_mags[i] - Ag
		rmag_dust = r_mags[i] - Ar
		imag_dust = i_mags[i] - Ai
		#colours
		g_i = gmag_dust - imag_dust
		g_r = gmag_dust - rmag_dust
		#k-corrected magnitudes
		lum_mpc = cosmo.luminosity_distance(z_clus[i]) #bcg z
		lum_mpc = lum_mpc.value
		Mg = gmag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
		Mr = rmag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
		Mi = imag_dust - 5*(np.log10(lum_mpc*(10**6))-1)
		kcor_r = calc_kcor('r',z_clus[i],'g - r',g_r)
		kcor_g = calc_kcor('g',z_clus[i],'g - r',g_r)
		kcor_i = calc_kcor('i',z_clus[i],'g - i',g_i)
		k_cor_Mr = Mr - kcor_r
		k_cor_Mi = Mi - kcor_i
		k_cor_Mg = Mg - kcor_g
		kcor_mr = rmag_dust - kcor_r
		kcor_mg = gmag_dust - kcor_g
		kcor_mi = imag_dust - kcor_i
		mass = 1.15 + 0.70*(kcor_mg -  kcor_mi) - 0.4*(k_cor_Mi) #restframe colours. Using Taylor et al., 2015 relation rather than 2011
		#Output
		Ag_out.append(Ag)
		Ar_out.append(Ar)
		Ai_out.append(Ai)
		k_cor_mr_out.append(kcor_mr)
		k_corr_mg_out.append(kcor_mg)
		k_corr_mr_out.append(kcor_mr)
		k_corr_mi_out.append(kcor_mi)
		k_corr_Mi_out.append(k_cor_Mi)
		k_corr_Mg_out.append(k_cor_Mg)
		k_corr_Mr_out.append(k_cor_Mr)
		mass_out.append(mass)
		kcor_r_out.append(kcor_r)
		kcor_i_out.append(kcor_i)
		kcor_g_out.append(kcor_g)
		ref_id.append(refname[i])
	else:
		Ag_out.append(np.nan)
		Ar_out.append(np.nan)
		Ai_out.append(np.nan)
		k_cor_mr_out.append(np.nan)
		k_corr_mg_out.append(np.nan)
		k_corr_mi_out.append(np.nan)
		k_corr_mr_out.append(np.nan)
		k_corr_Mi_out.append(np.nan)
		k_corr_Mg_out.append(np.nan)
		k_corr_Mr_out.append(np.nan)
		mass_out.append(np.nan)
		kcor_r_out.append(np.nan)
		kcor_i_out.append(np.nan)
		kcor_g_out.append(np.nan)
		ref_id.append(refname[i])

tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='25A', array=ref_id), \
	fits.Column(name='Ag', format='E', array=Ag_out), \
	fits.Column(name='Ar', format='E', array=Ar_out), \
	fits.Column(name='Ai', format='E', array=Ai_out), \
	fits.Column(name='kcor_g', format='E', array=kcor_g_out), \
	fits.Column(name='kcor_r', format='E', array=kcor_r_out), \
	fits.Column(name='kcor_i', format='E', array=kcor_i_out), \
	fits.Column(name='mg_kdust', format='E', array=k_corr_mg_out), \
	fits.Column(name='mr_kdust', format='E', array=k_corr_mr_out), \
	fits.Column(name='mi_kdust', format='E', array=k_corr_mi_out), \
	fits.Column(name='Mg_kdust', format='E', array=k_corr_Mg_out), \
	fits.Column(name='Mr_kdust', format='E', array=k_corr_Mr_out), \
	fits.Column(name='Mi_kdust', format='E', array=k_corr_Mi_out), \
	fits.Column(name='MStellar_Taylor_kdust', format='E', array=mass_out), \
										])
tbhdu.writeto(url+'kcorr_and_dustcorr_stellarmass.fits')
################################################################################
