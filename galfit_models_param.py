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
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import matplotlib.mlab as mlab
import math
from numpy import arccos
from numpy import cos
from numpy import sin
from astropy.cosmology import z_at_value
#from matplotlib.pyplot import cm
#import itertools
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# from fit_ellipse import fit_ellipse
from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic2D
import astropysics as pysics
from astropysics.models import SersicModel
import csv
from astropy.convolution import convolve_fft, convolve
from shutil import copyfile
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
################################################################################
url = '/home/kate/Desktop/'
folder_1 = 'clusters_for_models/'
folder_1a = folder_1+'folder_1/'
folder_2 = 'SIGMA/'
line_replace = [3,5,11,20,21,22,23,24,25]
gala_file = open(url+folder_1a+'galfit.feedme')
galafile = []
for line in gala_file.readlines():
	line = line.strip()
	galafile.append([line])
ref_strings = ['max','min','median']
bands = ['g','r','i']
# bins = fits.open(url+folder_1+'data_params_for_model_sample.fits')
# bins = bins[1].data

################################################################################
#get/convert psfs to counts
refnames = fits.open(url+folder_1+'model_refnames.fits')
refnames = refnames[1].data
ref = refnames.field('REFNAME')
# for i in range(0,len(bands)):
# 	for g in range(0,len(refnames)):
# 		psf = fits.open(url+folder_2+intstr(refnames[g])+'/repo/'+intstr(refnames[g])+'_sdss_'+bands[i]+'_psf.fits')
# 		clus_params = fits.open(url+'SIGMA/'+intstr(refnames[g])+'/sigfit/'+intstr(refnames[g])+'_sdss_'+bands[i]+'_modS.sigma.fits')
# 		clus_params_2 = fits.open(url+'SIGMA/'+intstr(refnames[g])+'_'+bands[i]+'.fits')
# 		counttonmgy = clus_params[0].header['NMGY']
# 		counttonmgy_2 = clus_params_2[0].header['NMGY']
# 		print counttonmgy_2,counttonmgy,bands[i],intstr(refnames[g])+'_sdss_'+bands[i]+'_modS.sigma.fits'
# 		psf = psf[0].data
# 		psf = psf/counttonmgy
# 		fits.writeto(url+folder_1a+intstr(refnames[g])+'_'+bands[i]+'_counts_psf.fits',np.array(psf))

#generate the feedme files
scale = cosmo.kpc_proper_per_arcmin(0.2) #shift to 0.2
scale = scale.value
Re = [10.0,20.0,30.0,50.0,100.0,150.0] #Re input (arcsec) Re_count,
Re_true = (np.array(Re)/scale)*(60.0/0.396)
Re_true = np.array(Re_true).tolist()
Re = Re_true
n_sersic = [1.0,2.0,4.0,8.0,12.0]
mag = [13.0,16.0,19.0]
ar = 0.66 #just keep this the same
pa = 50 #ditto pa
# median_name = ['74','182','112','163','65','345','73','56']
# x_0 = refnames.field('X')
# y_0 = refnames.field('Y')
# #Re n mag ar
# sigma_infile = []
# header = ['REFNAME','X','Y']
# for i in range(0,len(bands)):
# 	# sigma_infile = []
# 	for j in range(0,len(Re)):
# 		for k in range(0,len(n_sersic)):
# 			for l in range(0,len(mag)):
# 				for p in range(0,len(median_name)):
# 					for q in range(0,len(ref)):
# 						if np.str(ref[q]) == median_name[p]:
# 							clus_params = fits.open(url+'SIGMA/main_run_bcg_v2_sample/'+intstr(ref[q])+'/sigfit/'+intstr(ref[q])+'_sdss_'+bands[i]+'_modS.sigma.fits')
# 							zp = clus_params[0].header['MAGZP']
# 							input_3 = ['B) '+intstr(ref[q])+'_'+intstr(j)+'_'+intstr(k)+'_'+intstr(l)+'_'+bands[i]+'.fits           # Output data image block']
# 							input_5 = ['D) '+intstr(ref[q])+'_'+bands[i]+'_counts_psf.fits                #        # Input PSF image and (optional) diffusion kernel']
# 							input_11 = ['J) '+np.str(zp)+'              # Magnitude photometric zeropoint']
# 							input_20 = ['1) '+np.str(x_0[q])+' '+np.str(y_0[q])+'   1 1    # position x, y        [pixel]']
# 							input_21 = ['3) '+np.str(mag[l])+'        1       # total magnitude    ']
# 							input_22 = ['4) '+np.str(Re[j])+'         1       #     R_e              [Pixels]']
# 							input_23 = ['5) '+np.str(n_sersic[k])+'       1       # Sersic exponent (deVauc=4, expdisk=1)  ']
# 							input_24 = ['9) '+np.str(ar)+'       1       # axis ratio (b/a)   ']
# 							input_25 = ['10) '+np.str(pa)+'       1       # position angle (PA)  [Degrees: Up=0, Left=90]']
#
# 							replace(url+folder_1a+'galfit.feedme',url+folder_1a+'galfit_temp.feedme',galafile[3][0],input_3[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[5][0],input_5[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[11][0],input_11[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[20][0],input_20[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[21][0],input_21[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[22][0],input_22[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[23][0],input_23[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1a+'galfit_temp.feedme',galafile[24][0],input_24[0])
# 							replace(url+folder_1a+'galfit_temp.feedme',url+folder_1+'sigma_paramspace_infiles/galfit_'+intstr(ref[q])+'_'+intstr(j)+'_'+intstr(k)+'_'+intstr(l)+'_'+bands[i]+'.feedme',galafile[25][0],input_25[0])
# 							remove(url+folder_1a+'galfit_temp.feedme')
# 							sigma_infile.append([intstr(ref[q])+'_'+intstr(j)+'_'+intstr(k)+'_'+intstr(l)+'_'+bands[i],np.str(np.round(x_0[q],8)),np.str(np.round(y_0[q],8)),Re[j],n_sersic[k],mag[l]])
# np.savetxt(url+folder_1+'sigma_paramspace_infiles/sigma_infile_allparameters.csv', sigma_infile, delimiter=',',fmt='%s')   # X is an array
# print 'Infiles created! Models being made....'
# # ##############################################################################
# #Run galfit to make models
# refnames = np.loadtxt(url+folder_1+'sigma_paramspace_infiles/sigma_infile_allparameters.csv',comments="#", delimiter=",", unpack=False,dtype='str')
# refnames = refnames[:,0]
# # print refnames[0]
# print intstr(len(refnames))+' to be modelled.'
# t = time.time()
# for q in range(0,len(refnames)):
# 	filename = 'galfit_'+refnames[q]+'.feedme'
# 	doFitting(filename)
# 	print time.time()-t
# print 'Models made, adding to fields...'
# ################################################################################
#add models to original images
orig_image_url = '/Desktop/SIGMA/main_run_bcg_v2_sample/'
median_name = ['74','182','112','163','65','345','73','56']
n_sersic = [1.0,2.0,4.0,8.0,12.0]
mag = [13.0,16.0,19.0]
ar = 0.66 #just keep this the same
pa = 50 #ditto pa
for g in range(0,len(median_name)): # generates model profile parameters
	for q in range(0,len(bands)): #for all bands, apply precisely the same profile (bar the magnitude)
		clus_params = fits.open(url+'SIGMA/main_run_bcg_v2_sample/'+intstr(median_name[g])+'/sigfit/'+intstr(median_name[g])+'_sdss_'+bands[q]+'_modS.sigma.fits')
		prihdr= clus_params[0].header
		zp = clus_params[0].header['MAGZP']
		counttonmgy = clus_params[0].header['NMGY']
		cluster = fits.open(url+'SIGMA/main_run_bcg_v2_sample/'+intstr(median_name[g])+'_'+bands[q]+'.fits',ignore_missing_end=True)
		cluster = cluster[0].data
		for j in range(0,len(Re)):
			for k in range(0,len(n_sersic)):
				for l in range(0,len(mag)):
					model = fits.open(url+folder_1+'sigma_paramspace_infiles/'+intstr(median_name[g])+'_'+intstr(j)+'_'+intstr(k)+'_'+intstr(l)+'_'+bands[q]+'.fits')
					model = model[0].data
					combined = (cluster + (model*counttonmgy))/counttonmgy #spits out model and cluster in ADU
					fits.writeto(url+folder_1+'sigma_paramspace_infiles/'+intstr(median_name[g])+'_'+intstr(j)+'_'+intstr(k)+'_'+intstr(l)+'_'+bands[q]+'_combined_galfit.fits',np.array(combined),header=prihdr)
# ################################################################################
