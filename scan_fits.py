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
#import itertools
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

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

def join_table(table1,table2): #note at least one needs object id
	catalog = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
	match = []
	for i in range(0,len(file_id)):
		ra1 = RA[i]
		dec1 = DEC[i]
		c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
		idx1, d2d1, d3d1 = c.match_to_catalog_sky(catalog_xclass)
		if d2d1.degree < 0.5:
			match.append([file_id[i],RA[i],DEC[i],idx1,np.float(d2d1.degree)])

def extract_fits_table(file_location):
	datafile = fits.open(file_location)
	datafile = datafile[1].data
	return datafile ;

def sdss_celestial_to_x_y(imfile, ra, dec):
	hdu = fits.open(imfile)
	hdu[0].header.rename_keyword('RADECSYS', 'RADESYS') #This is an annoying error from the sdss files that fucked you over!
	w = WCS(hdu[0].header)
	pix = w.wcs_world2pix(ra,dec,1) #YES. YES. YEEEEEEESSSS!!!!
	return pix;

def fits_save(headers,table,fileloc): #this is very simplistic
	colarray = []
	for i in range(0,len(headers)):
		data = table[i]
		header= headers[i]
		if type(data[0]) == 'float':
			data_form = 'E'
		elif type(data[0]) == 'string':
			data_form = '20A'
		col = fits.Column(name=header, format=data_form, array=data)
		colarray.append(col)
	tbhdu = fits.BinTableHDU.from_columns(colarray)
	tbhdu.writeto(fileloc)
	return ;

def gapper_biweight(v): #returns nans
	limit = 10
	n_obj = len(v)
	ksp = np.nan
	v_sorted = []
	sigma_cut = []
	gapper = np.nan
	biweight = np.nan
	obj_flag = np.nan
	if n_obj == 1: #if there's one spec_z, chuck it out
		flag = 0 #0 = do not use
		n_obj = 1
	elif n_obj == 0:
		flag = 0
		n_obj = 0
	else: #if not, then:
		v_median = np.median(v) #take the median
		v_cut = []
		g = 0
		for i in range(0,len(v)): #for the remaining objects
			if v[i] < v_median + 3000 and v[i] > v_median - 3000: #clip at +/- 5000kms about median
				v_cut.append(v[i])
				g = g + 1 #keep count
		if g < limit: #if there's less than 10 objects left, chuck
			flag = 0
			n_obj = len(v_cut)
		else: #sigma clip up to convergence
			v_3 = v_cut
			sigma_cut = asts.sigma_clip(v_3,sigma=3,iters=None)
			sigma_cut = MaskedArray.compressed(sigma_cut)
			#print len(v_3), len(sigma_cut)
		if len(sigma_cut) < limit:
			gapper = np.nan
			biweight = np.nan
			flag = 0
			n_obj = len(sigma_cut)
		else:
			#GAPPER (CORRECT THIS TIME, WE HOPE)
			n_obj = len(sigma_cut)
			v_sorted = sorted(sigma_cut,reverse=True) #sort the remaining values
			#print v_sorted
			product = []
			coeff = (np.sqrt(np.pi))/(n_obj*(n_obj - 1.0))  #coefficient in estimator
			for i in range(0,len(v_sorted)-1):
				j = i+1
				gap = v_sorted[i] - v_sorted[j]
				weight = i*(n_obj - i)
				product.append(weight*gap)
			gapper = coeff*(sum(product))
			#BIWEIGHT (FOR COMPARISON)
			biweight = asts.biweight_midvariance(sigma_cut,M=np.median(sigma_cut))
			#print biweight,  gapper
			flag = 1
			# kstest_v = kstest(v_sorted,stats.norm.cdf, args=(np.mean(v_sorted),np.std(v_sorted))) #compares v with a normal distribution
			# ksp = kstest_v[1]
			if n_obj < 15:
				obj_flag = 'GAPPER'
			else:
				obj_flag = 'BIWEIGHT'
	return gapper, biweight, flag, n_obj, v_sorted, obj_flag;

def gapper_biweight_err(v,flag,n,biweight_true,gapper_true): #returns nans
	limit = 10
	n_obj = len(v)
	upper_biweight = np.nan
	lower_biweight=np.nan
	upper_gapper=np.nan
	lower_gapper=np.nan
	gapper_sigma = np.nan
	biweight_sigma=np.nan
	gapper_errfun = []
	biweight_errfun = []
	biweight = np.nan
	gapper=np.nan
	#print v
	#print flag
	if flag == 1: #if not, then
		#print biweight
		boot_resamp = asts.bootstrap(np.array(v),n)
		for p in range(0,n-1):
			dist = boot_resamp[p,:]
			#print dist - v
			#print dist - v
			#time.sleep(5)
			dist_median = np.median(dist)
			samp_cut = []
			#samp_cut = dist
			g = 0
			for c in range(0,len(dist)): #for the remaining objects
				if dist[c] < dist_median + 3000 and dist[c] > dist_median - 3000: #clip at +/- 3000kms about median
					samp_cut.append(dist[c])
					g = g + 1 #keep count
			#print g
			if g < limit: #if there's less than 10 objects left, chuck
				gapper = np.nan
				biweight = np.nan
				samp_cut = []
			else: #sigma clip up to convergence
				v_3 = samp_cut

			if len(samp_cut) < limit:
				sigma_cut = []
				biweight = np.nan
				gapper = np.nan
			elif len(samp_cut) >= limit:
				sigma_cut = asts.sigma_clip(v_3,sigma=3,iters=None)
				sigma_cut = MaskedArray.compressed(sigma_cut)
				#time.sleep(5)
				#print len(sigma_cut),len(v)
			if len(sigma_cut) < limit:
				biweight = np.nan
				gapper=np.nan
			else:
				biweight = asts.biweight_midvariance(sigma_cut,M=np.median(sigma_cut))
				n_obj = len(sigma_cut)
				v_sorted = sorted(sigma_cut,reverse=True) #sort the remaining values
				product = []
				coeff = (np.sqrt(np.pi))/(n_obj*(n_obj - 1.0))  #coefficient in estimator
				for i in range(0,len(v_sorted)-1):
					j = i+1
					gap = v_sorted[i] - v_sorted[j]
					weight = i*(n_obj - i)
					product.append(weight*gap)
					gapper = coeff*(sum(product))
			if biweight != np.nan:
				biweight_errfun.append(biweight)
			if gapper != np.nan:
				gapper_errfun.append(gapper)
		biweight_errfun = np.array(biweight_errfun)
		gapper_errfun = np.array(gapper_errfun)
		#print len(biweight_errfun),len(gapper_errfun)
		biweight_errfun = biweight_errfun[~np.isnan(biweight_errfun)]
		gapper_errfun = gapper_errfun[~np.isnan(gapper_errfun)]
		#print len(biweight_errfun),len(gapper_errfun)
		#gapper_err = stats.norm.interval(0.68, loc=np.mean(gapper_errfun), scale=np.std(gapper_errfun)/np.sqrt(len(gapper_errfun))) #assume gaussian
		biweight_conf = stats.norm.interval(0.68, loc=biweight_true, scale=asts.mad_std(biweight_errfun)/np.sqrt(len(biweight_errfun))) #assume gaussian
		#biweight_sigma=asts.mad_std(biweight_errfun)
		upper_biweight = biweight_conf[1] - biweight_true
		lower_biweight = biweight_true - biweight_conf[0]

		gapper_conf = stats.norm.interval(0.68, loc=gapper_true, scale=asts.mad_std(gapper_errfun)/np.sqrt(len(gapper_errfun))) #assume gaussian

		upper_gapper = gapper_conf[1] - gapper_true
		lower_gapper = gapper_true - gapper_conf[0]
		biweight_sigma=asts.mad_std(biweight_errfun)
		gapper_sigma=asts.mad_std(gapper_errfun)
		# print biweight_conf,gapper_conf,biweight_true,np.median(biweight_errfun)
		# #print upper_biweight,lower_biweight,upper_gapper,lower_gapper
	return upper_biweight,lower_biweight,upper_gapper,lower_gapper,biweight_sigma,gapper_sigma;

def biweight_z(v,z,z_phot): #returns nans
	limit = 10
	n_obj = len(v)
	ksp = np.nan
	v_sorted = []
	index_obj = np.nan
	biweight = np.nan
	sigma_cut = []
	if n_obj == 1: #if there's one spec_z, chuck it out
		flag = 0 #0 = do not use
		n_obj = 1
	elif n_obj == 0:
		flag = 0
		n_obj = 0
	else: #if not, then:
		v_median = np.median(v) #take the median
		v_cut = []
		g = 0
		for i in range(0,len(v)): #for the remaining objects
			if v[i] < v_median + 3000 and v[i] > v_median - 3000: #clip at +/- 5000kms about median
				v_cut.append(z[i])
				g = g + 1 #keep count
		if g < limit: #if there's less than 10 objects left, chuck
			flag = 0
			n_obj = len(v_cut)
		else: #sigma clip up to convergence
			v_3 = v_cut
			sigma_cut = asts.sigma_clip(v_3,sigma=3,iters=None)
			sigma_cut = MaskedArray.compressed(sigma_cut)
			# sigma_mask = [z == sigma_cut]
			# print sigma_mask
			# time.sleep(5)
		if len(sigma_cut) < limit:
			biweight = np.nan
			flag = 0
			n_obj = len(sigma_cut)
		else:
			#GAPPER (CORRECT THIS TIME, WE HOPE)
			n_obj = len(sigma_cut)
			#BIWEIGHT (FOR COMPARISON)
			biweight = asts.biweight_location(sigma_cut,M=np.median(sigma_cut))
			flag = 1
			# index_obj = np.array(indexval)[sigma_mask]
			# index_obj = np.array(index_obj).tolist()
			#print index_obj
			# time.sleep(5)
			# print biweight, z_phot
	return biweight, sigma_cut,	flag;

def biweight_err_z(z,n,flag,z_true): #returns nans
	limit = 10
	H_0 = 70
	n_obj = len(z)
	upper_biweight = np.nan
	lower_biweight=np.nan
	biweight_sigma=np.nan
	biweight_errfun = []
	biweight = np.nan
	#print v
	#print flag
	if flag ==1:
		boot_resamp = asts.bootstrap(np.array(z),n)
		for p in range(0,n-1):
			dist = boot_resamp[p,:]
			#print dist
			#print dist - v, len(dist),len(v)
			#time.sleep(5)
			dist_median = np.median(dist)
			comoving = cosmo.comoving_distance(dist_median)
			comoving_median = (comoving.value)*H_0
			samp_cut = []
			#samp_cut = dist
			g = 0
			for c in range(0,len(dist)): #for the remaining objects
				comoving = cosmo.comoving_distance(dist[c])
				comoving = (comoving.value)*H_0
				if comoving < comoving_median + 3000 and comoving > comoving_median - 3000: #clip at +/- 5000kms about median
					samp_cut.append(dist[c])
					g = g + 1 #keep count
					#print g
			if g < limit: #if there's less than 10 objects left, chuck
				biweight = np.nan
				samp_cut = []
			else: #sigma clip up to convergence
				v_3 = samp_cut
			if len(samp_cut) < limit:
				sigma_cut = []
				biweight = np.nan
			elif len(samp_cut) >= limit:
				sigma_cut = asts.sigma_clip(v_3,sigma=3,iters=None)
				sigma_cut = MaskedArray.compressed(sigma_cut)
			#time.sleep(5)
			#print len(sigma_cut),len(v)
			if len(sigma_cut) < limit:
				biweight = np.nan
			else:
				biweight = asts.biweight_location(sigma_cut,M=np.median(sigma_cut))
			if biweight != np.nan:
				biweight_errfun.append(biweight)
		biweight_errfun = np.array(biweight_errfun)
		biweight_errfun = biweight_errfun[~np.isnan(biweight_errfun)]
		#gapper_err = stats.norm.interval(0.68, loc=np.mean(gapper_errfun), scale=np.std(gapper_errfun)/np.sqrt(len(gapper_errfun))) #assume gaussian
		biweight_conf = stats.norm.interval(0.68, loc=z_true, scale=asts.mad_std(biweight_errfun)/np.sqrt(len(biweight_errfun))) #assume gaussian
		biweight_sigma=asts.mad_std(biweight_errfun)
		upper_biweight = biweight_conf[1] - z_true
		lower_biweight = z_true - biweight_conf[0]
		biweight_sigma=asts.mad_std(biweight_errfun)
		#print upper_biweight,lower_biweight,np.median(biweight_errfun),z_true
		# time.sleep(2)
		# print sigma_cut
		# time.sleep(2)
	return upper_biweight,lower_biweight,biweight_sigma;

################################################################################
bands = ['g','r','i']

url = '/home/kate/Desktop/'

sample_old = extract_fits_table(url+'New_BCGs/image_data/confirmed/orig_sample_first_run.fits')
sample_new = extract_fits_table(url+'New_BCGs/image_data/confirmed/sample_first_fit_run.fits')

f_id_old = sample_old.field('CLUS_ID_1a')
f_id_new = sample_new.field('CLUS_REF')

spurious_fits_new = []
for i in range(0,len(f_id_new)):
	z = sample_new.field('Z_BCG')[i]
	scale = cosmo.kpc_proper_per_arcmin(z)
	scale = scale.value
	band_rej = [f_id_new[i]]
	for g in range(0,len(bands)):
		re = sample_new.field('GAL_sdss_'+bands[g]+'_modS_C1_RE')[i]
		n = sample_new.field('GAL_sdss_'+bands[g]+'_modS_C1_RE')[i]
		re_kpc = scale*(re*(0.396/60))
		if re_kpc > 100 or n < 1 or n > 14:
			band_rej.extend([bands[g]])
		else:
			band_rej.extend(['n'])
	spurious_fits_new.append(band_rej)

g_fit = []
r_fit = []
i_fit = []
for i in range(0,len(spurious_fits_new)):
	g_fit.append(spurious_fits_new[i][1])
	r_fit.append(spurious_fits_new[i][2])
	i_fit.append(spurious_fits_new[i][3])
tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=f_id_new), \
		fits.Column(name='refit_g', format='20A', array=g_fit), \
		fits.Column(name='refit_r', format='20A', array=r_fit), \
		fits.Column(name='refit_i', format='20A', array=i_fit)])

if  os.path.isfile(url+'New_BCGs/image_data/confirmed/refits_new_sample.fits') == True:
	os.remove(url+'New_BCGs/image_data/confirmed/refits_new_sample.fits')
tbhdu.writeto(url+'New_BCGs/image_data/confirmed/refits_new_sample.fits')

spurious_fits_old = []
for i in range(0,len(f_id_old)):
	z = sample_old.field('Z_ALL')[i]
	scale = cosmo.kpc_proper_per_arcmin(z)
	scale = scale.value
	band_rej = [f_id_old[i]]
	for g in range(0,len(bands)):
		re = sample_old.field('GAL_sdss_'+bands[g]+'_modS_C1_RE')[i]
		n = sample_old.field('GAL_sdss_'+bands[g]+'_modS_C1_RE')[i]
		re_kpc = scale*(re*(0.396/60))
		if re_kpc > 100 or n < 1 or n > 14:
			band_rej.extend([bands[g]])
		else:
			band_rej.extend(['n'])
	spurious_fits_old.append(band_rej)

g_fit = []
r_fit = []
i_fit = []
# f_id_old = [float(i) for i in f_id_old]
for i in range(0,len(spurious_fits_old)):
	g_fit.append(spurious_fits_old[i][1])
	r_fit.append(spurious_fits_old[i][2])
	i_fit.append(spurious_fits_old[i][3])
tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=f_id_old), \
		fits.Column(name='refit_g', format='20A', array=g_fit), \
		fits.Column(name='refit_r', format='20A', array=r_fit), \
		fits.Column(name='refit_i', format='20A', array=i_fit)])

if os.path.isfile(url+'SIGMA/main_run_bcg_v2_sample/refits_old_sample.fits') == True:
	os.remove(url+'SIGMA/main_run_bcg_v2_sample/refits_old_sample.fits')
tbhdu.writeto(url+'SIGMA/main_run_bcg_v2_sample/refits_old_sample.fits')
################################################################################
#visinspec
image_url_old = url+'SIGMA/main_run_bcg_v2_sample/'
image_url_new = url+'/New_BCGs/image_data/confirmed/SIGMA/'
