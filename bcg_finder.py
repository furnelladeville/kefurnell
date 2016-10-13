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
	return match ;

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
			data_form = '25A'
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
#Testing file saving method for neighbour z-shifts
url = '/home/kate/Desktop/New_BCGs/'
folder = 'image_data/'
#v_disps = extract_fits_table(url+'confirmed_candidates.fits')
v_disps = extract_fits_table('/home/kate/Desktop/SIGMA/main_run_bcg_v2_sample/final_fitted_sample_10_z.fits')
clus_id = v_disps.field('CLUS_ID_1a')
z_shifts = v_disps.field('ALLZ_NOQSO')
n_z = v_disps.field('NHASZ')
r200c = v_disps.field('R200c/deg')
ra_obj = v_disps.field('ALLRA')
dec_obj = v_disps.field('ALLDEC')
ra_bcg = v_disps.field('RA_BCG') #Use the BCG for now. It's the best you can really do. Ask Chris?
dec_bcg = v_disps.field('DEC_BCG')
z_phot= v_disps.field('Z_LAMBDA') #photo-z, going to have to use this to make comparison
# bcg_save = []
# for i in range(0, len(clus_id)):
# 	z = np.array(z_shifts[i])
# 	z = z[0:n_z[i]]
# 	z = np.array(z).tolist()
# 	diff = -np.ones(np.int(60-n_z[i]),dtype=np.float)
# 	diff = np.array(diff).tolist()
# 	z = z + diff
# 	bcg_save.append(z)
# tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=clus_id), \
# 	fits.Column(name='Z_BCG', format='60E', array=bcg_save)])
# tbhdu.writeto(url+'test.fits')
################################################################################
# H_0 = 70
# c2 = 3*10**5
# Omega_l = 0.7
# Omega_m = 0.3
# z_bwt = []
# z_bwt_err = []
# sigma_save = []
# sigma_err_save = []
# gals_kept = []
# m200_save = []
# m200_err_save =[]
# r200_save = []
# r200_err_save = []
# sigma_flag = []
# n_kept = []
# start_time = time.time()
# for i in range(0,len(clus_id)):
# 	ra_clus = ra_bcg[i]
# 	dec_clus = dec_bcg[i]
# 	objs = z_shifts[i][0:n_z[i]] #store the indices to extract out the values you keep later
# 	ra = ra_obj[i][0:n_z[i]]
# 	dec = dec_obj[i][0:n_z[i]]
# 	z_obj = []
# 	index_obj = np.array(range(0,n_z[i]))
# 	v_objs_comoving = []
# 	#print index
# 	for q in range(0,len(objs)):
# 		sep = np.arccos(cos(90-dec[q])*cos(90-dec_clus) + sin(90-dec[q])*sin(90-dec_clus)*cos(ra[q]-ra_clus))
# 		if sep < r200c[i]:
# 			comoving = cosmo.comoving_distance(objs[q])
# 			comoving = comoving.value
# 			v_objs_comoving.append(H_0*comoving)
# 			z_obj.append(objs[q])
# 	zbwt = biweight_z(v_objs_comoving,z_obj,z_phot[i])
# 	zbwt_err = biweight_err_z(z_obj,10000,zbwt[2],zbwt[0])
# 	#biweight, sigma_cut, flag, index_obj;
# 	v2 = ((np.array(zbwt[1]) - zbwt[0])/(1+zbwt[0]))*c2
# 	sig_clus = gapper_biweight(v2)
# 	sig_clus_err = gapper_biweight_err(v=v2,flag=sig_clus[2],n=10000, biweight_true=sig_clus[1],gapper_true=sig_clus[0])
# 	if sig_clus[2] == 1 and sig_clus[5] == 'GAPPER': #Taken from Rose Finn's 2005 paper - she recommends using X-ray gas mass as a proxy if cluster is not relaxed
# 		r200 = (1.73)*((sig_clus[0])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		r200_err = (1.73)*((sig_clus_err[5])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		m200 = (1.0/0.7)*(1.2*10**15.0)*((sig_clus[0]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		m200_err = (1.0/0.7)*(1.2*10**15.0)*((sig_clus_err[5]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		sigma = sig_clus[0]
# 		sigma_err = sig_clus_err[5]
# 		used_mask = np.in1d(np.array(objs),np.array(zbwt[1]))
# 		used_obj = index_obj[used_mask]
# 		numberkept = len(used_obj)
# 		used_obj = np.array(used_obj).tolist() + np.array(-np.ones(60-len(used_obj))).tolist()
# 		best = sig_clus[5]
# 	elif sig_clus[2] == 1 and sig_clus[5] == 'BIWEIGHT': #Taken from Rose Finn's 2005 paper - she recommends using X-ray gas mass as a proxy if cluster is not relaxed
# 		r200 = (1.73)*((sig_clus[1])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		r200_err = (1.73)*((sig_clus_err[4])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		m200 = (1.0/0.7)*(1.2*10**15.0)*((sig_clus[1]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		m200_err = (1.0/0.7)*(1.2*10**15.0)*((sig_clus_err[4]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
# 		sigma = sig_clus[1]
# 		sigma_err = sig_clus_err[4]
# 		used_mask = np.in1d(np.array(objs),np.array(zbwt[1]))
# 		used_obj = index_obj[used_mask]
# 		numberkept = len(used_obj)
# 		used_obj = np.array(used_obj).tolist() + np.array(-np.ones(60-len(used_obj))).tolist()
# 		best = sig_clus[5]
# 	else:
# 		r200 = np.nan
# 		r200_err = np.nan
# 		m200 = np.nan
# 		m200_err=np.nan
# 		sigma = np.nan
# 		sigma_err = np.nan
# 		numberkept = np.nan
# 		used_obj = np.array(-np.ones(60)).tolist()
# 		best = 'NONE'
# 	gals_kept.append(used_obj)
# 	z_bwt.append(zbwt[0])
# 	z_bwt_err.append(zbwt_err[2])
# 	sigma_save.append(sigma)
# 	sigma_err_save.append(sigma_err)
# 	m200_save.append(m200)
# 	m200_err_save.append(m200_err)
# 	r200_save.append(r200)
# 	r200_err_save.append(r200_err)
# 	sigma_flag.append(best)
# 	n_kept.append(numberkept)
# 	print r200,m200,zbwt[0],sigma,sigma_err,used_obj,best,numberkept
# 	print("--- %s seconds ---" % (time.time() - start_time))
# tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=clus_id), \
# 	fits.Column(name='SIGMA_BEST', format='E', array=sigma_save), \
# 	fits.Column(name='SIGMA_BEST_ERR', format='E', array=sigma_err_save), \
# 	fits.Column(name='Z_BWT', format='E', array=z_bwt), \
# 	fits.Column(name='Z_BWT_ERR', format='E', array=z_bwt_err), \
# 	fits.Column(name='M200_BEST', format='E', array=m200_save), \
# 	fits.Column(name='M200_BEST_ERR', format='E', array=m200_err_save), \
# 	fits.Column(name='R200_BEST', format='E', array=r200_save), \
# 	fits.Column(name='R200_BEST_ERR', format='E', array=r200_err_save), \
# 	fits.Column(name='N_KEPT', format='E', array=n_kept), \
# 	fits.Column(name='BEST_SIGMA', format='20A', array=sigma_flag), \
# 	fits.Column(name='KEPT_INDEX', format='60E', array=gals_kept)])
# # tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=clus_id), \
# # 	fits.Column(name='Z_BWT_ERR', format='E', array=z_bwt_err)])
# tbhdu.writeto('/home/kate/Desktop/SIGMA/main_run_bcg_v2_sample/v_disp_recalc_fullfitted.fits')
################################################################################
#test sdss celestial
# n_missing = 0
url = '/home/kate/Desktop/New_BCGs/'
folder = 'image_data/'
# imcoords = extract_fits_table(url+'image_data/new_sample_for_inspection.fits')
# refname2 = imcoords.field('CLUS_ID')
# alt_xcoords = extract_fits_table(url+folder+'alt_opt/alt_optcoords_check.fits')
# filename = alt_xcoords.field('filename')
# refnames = alt_xcoords.field('CLUS_REF')
# x_x = []
# y_x = []
# ra = []
# dec = []
# x_bcg = []
# y_bcg = []
# missing_flag = []
# clus_idlist = []
# z_bcg = []
# for i in range(0,len(filename)):
# 	filemask = [refnames == refnames[i]]
# 	#print filemask
# 	file_id = filename[filemask]
# 	#print file_id
# 	filemask = [refname2 == refnames[i]] #GET STUFF OUT OF MAIN FOLDER
# 	#print filemask
# 	refname = refname2[filemask]
# 	CODEX_RA = imcoords.field('RA_2')[filemask]
# 	CODEX_DEC = imcoords.field('DEC_2')[filemask]
# 	RA_BCG = imcoords.field('ALLRA')[filemask]
# 	DEC_BCG = imcoords.field('ALLDEC')[filemask]
# 	FLAG_BCG = imcoords.field('ALLBCGFLAG')[filemask]
# 	items = len(RA_BCG) #doesn't matter too much
# 	n_z = imcoords.field('NHASZ')[filemask]
# 	Z_BCG = imcoords.field('ALLZ_NOQSO')[filemask]
# 	conv_codex = sdss_celestial_to_x_y(url+folder+'alt_opt/'+file_id[0], ra=CODEX_RA, dec=CODEX_DEC)
# 	if 1.0 in FLAG_BCG:
# 		flagmask = [FLAG_BCG == 1]
# 		RA_BCG1 = RA_BCG[flagmask]
# 		DEC_BCG1 = DEC_BCG[flagmask]
# 		Z_BCG1 = Z_BCG[flagmask]
# 		conv = sdss_celestial_to_x_y(url+folder+'alt_opt/'+file_id[0], ra=RA_BCG1, dec=DEC_BCG1)
# 		xval = conv[0]
# 		yval = conv[1]
# 		missing = 'n'
# 	else:
# 		RA_BCG1 = np.nan
# 		DEC_BCG1 = np.nan
# 		Z_BCG1 = np.nan
# 		missing = 'y'
# 		xval = np.nan
# 		yval = np.nan
# 		n_missing = n_missing + 1
# 	#print file_id, RA_BCG1,DEC_BCG1
# 	#pixel_coords.append([refname,RA_BCG1,DEC_BCG1,np.array(conv[0]).tolist(),np.array(conv[1]).tolist(),missing,np.array(conv_codex[0]).tolist(),np.array(conv_codex).tolist()])
# 	x_x.append(conv_codex[0])
# 	y_x.append(conv_codex[1])
# 	x_bcg.append(xval)
# 	y_bcg.append(yval)
# 	ra.append(RA_BCG1)
# 	dec.append(DEC_BCG1)
# 	missing_flag.append(missing)
# 	clus_idlist.append(refname)
# 	z_bcg.append(Z_BCG1)
# 	# bcg_data = []
# 	# for l in range(0,items):
# 	# 	if l < NHASZ and RA_BCG[l] >-500:
# 	# 		conv = sdss_celestial_to_x_y(url+folder+file_id, ra=RA_BCG[l], dec=DEC_BCG[l])
# 	# 		bcg_data.append[conv[0],conv[1]]
# #print pixel_coords, n_missing
# #print len(x_x), len(y_x), len(x_bcg)
# tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=clus_idlist), \
# 	fits.Column(name='MISSING', format='20A', array=missing_flag), \
# 	fits.Column(name='Z_BCG', format='E', array=z_bcg), \
# 	fits.Column(name='RA_BCG', format='E', array=ra), \
# 	fits.Column(name='DEC_BCG', format='E', array=dec), \
# 	fits.Column(name='X_BCG', format='E', array=x_bcg), \
# 	fits.Column(name='Y_BCG', format='E', array=y_bcg), \
# 	fits.Column(name='X_CODEX', format='E', array=np.array(x_x).tolist()), \
# 	fits.Column(name='Y_CODEX', format='E', array=np.array(y_x).tolist())])
# tbhdu.writeto(url+folder+'alt_opt/clus_alt_opterr.fits')
################################################################################
#PLOT ON JPEG IMAGES TO LOOK THROUGH
# bcg_data = extract_fits_table(url+folder+'alt_opt/clus_alt_opterr.fits')
# bcg_data_2 = extract_fits_table(url+folder+'new_bcg_data_xy.fits')
# all_data = extract_fits_table(url+folder+'new_sample_for_inspection.fits')
# clus_id2 = all_data.field('CLUS_ID')
# clus_id = bcg_data.field('CLUS_REF')
# x_codex = bcg_data.field('X_CODEX')
# y_codex = bcg_data.field('Y_CODEX')
# x_bcg = bcg_data.field('X_BCG')
# y_bcg = bcg_data.field('Y_BCG')
# z_bcg = bcg_data.field('Z_BCG')
# missing = bcg_data.field('MISSING')
# for i in range(0,len(clus_id)):
# 	if missing[i] == 'n':
# 		mask = [clus_id2 == clus_id[i]]
# 		z_all = all_data.field('ALLZ_NOQSO')[mask]
# 		z_all = z_all[0]
# 		#print z_all
# 		nmem = all_data.field('NHASZ')[mask]
# 		#print nmem
# 		x_bcg2 = bcg_data_2.field('X_BCG')[mask]
# 		y_bcg2 = bcg_data_2.field('Y_BCG')[mask]
# 		x_codex2 = bcg_data_2.field('X_CODEX')[mask]
# 		y_codex2 = bcg_data_2.field('y_CODEX')[mask]
#
# 		median_z = np.median(np.array(z_all[0:nmem]))
# 		#print median_z
# 		image = mpimage.imread(url+folder+'alt_opt/'+clus_id[i]+'.jpeg')
# 		image2 = mpimage.imread(url+folder+clus_id[i]+'.jpeg')
# 		plt.subplot(121)
# 		plt.imshow(image)
# 		plt.plot(x_bcg[i],1489-y_bcg[i], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='b',markeredgewidth=1.2)
# 		plt.plot(x_codex[i],1489-y_codex[i], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='r',markeredgewidth=1.2)
# 		plt.axis('tight')
# 		plt.subplot(122)
# 		plt.imshow(image2)
# 		plt.plot(x_bcg2,1489-y_bcg2, marker='s',markersize=8,markerfacecolor='none',markeredgecolor='b',markeredgewidth=1.2)
# 		plt.plot(x_codex2,1489-y_codex2, marker='s',markersize=8,markerfacecolor='none',markeredgecolor='r',markeredgewidth=1.2)
# 		plt.axis('tight')
# 		plt.show()
# 		print clus_id[i], z_bcg[i], median_z
# 	else:
# 		print clus_id[i]
################################################################################
#rename files
# filenames_renamed_already = ['1_6107','2_19458','2_19623','1_4021','1_507','1_6181','1_3024','1_1275']
# bcg_data = extract_fits_table(url+folder+'new_sample_for_inspection.fits')
# refname = bcg_data.field('CLUS_ID')
# run = bcg_data.field('RUN')
# camcol = bcg_data.field('CAMCOL')
# fields = bcg_data.field('FIELD')
# for i in range(0,len(refname)):
# 	if refname[i] not in filenames_renamed_already:
# 		#print fields[i]
# 		if np.int(fields[i]) <10:
# 			field = '000'+intstr(fields[i])
# 		elif np.int(fields[i]) >=10 and np.int(fields[i]) < 100:
# 			field = '00'+intstr(fields[i])
# 		elif np.int(fields[i]) >100:
# 			field = '0'+intstr(fields[i])
# 		r_id = 'frame-r-00'+intstr(run[i])+'-'+intstr(camcol[i])+'-'+field+'.fits'
# 		g_id = 'g/frame-g-00'+intstr(run[i])+'-'+intstr(camcol[i])+'-'+field+'.fits'
# 		i_id = 'i/frame-i-00'+intstr(run[i])+'-'+intstr(camcol[i])+'-'+field+'.fits'
# 		#print url+folder+r_id,refname[i]+'_r.fits'
# 		#time.sleep(5)
# 		print refname[i]
# 		# if os.path.isfile(url+folder+r_id) == True: #Checks to see if file exists in path. You don't even need the array at the beginning of the loop....
# 		# 	os.rename(url+folder+r_id,url+folder+refname[i]+'_r.fits')
# 		if os.path.isfile(url+folder+g_id) == True:
# 			os.rename(url+folder+g_id,url+folder+'g/'+refname[i]+'_g.fits')
# 		# if os.path.isfile(url+folder+i_id) == True:
# 		# 	os.rename(url+folder+i_id,url+folder+'i/'+refname[i]+'_i.fits')
################################################################################
#Take out coordinates from main file for ease
# full_sample = extract_fits_table(url+'new_sample_for_inspection.fits')
# samp_ids = full_sample.field('CLUS_ID')
# zall = full_sample.field('ALLZ_NOQSO')
# zallerr = full_sample.field('ALLZ_ERR_NOQSO')
# bcgflag = full_sample.field('ALLBCGFLAG')
# rabcg = full_sample.field('ALLRA')
# decbcg = full_sample.field('ALLDEC')
# confirmed = extract_fits_table(url+folder+'confirmed/confirmed_ids.fits')
# confirmed_id = confirmed.field('CLUS_ID') #these are just the cluster ids
# confirmed_rank = confirmed.field('BCGFLAG')
# RA_BCG = []
# DEC_BCG = []
# Z_BCG = []
# Z_ERR_BCG = []
# ID = []
# for i in range(0,len(confirmed_id)):
# 	clus_mask = [samp_ids == confirmed_id[i]]
# 	bcg_mask = [bcgflag[clus_mask] == confirmed_rank[i]]
# 	z = zall[clus_mask][bcg_mask]
# 	zerr = zallerr[clus_mask][bcg_mask]
# 	ra = rabcg[clus_mask][bcg_mask]
# 	dec = decbcg[clus_mask][bcg_mask]
# 	Z_BCG.append(z)
# 	Z_ERR_BCG.append(zerr)
# 	ID.append(confirmed_id[i])
# 	RA_BCG.append(ra)
# 	DEC_BCG.append(dec)
# tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='CLUS_REF', format='20A', array=ID), \
# 	fits.Column(name='Z_BCG', format='E', array=Z_BCG), \
# 	fits.Column(name='Z_ERR', format='E', array=Z_ERR_BCG), \
# 	fits.Column(name='RA_BCG', format='E', array=RA_BCG), \
# 	fits.Column(name='DEC_BCG', format='E', array=DEC_BCG)])
# tbhdu.writeto(url+folder+'confirmed/candidate_coords_redshifts.fits')
################################################################################
#move confirmed files around
#url = '/home/kate/Desktop/New_BCGs/'
url = '/home/kate/Desktop/'
folder = 'SIGMA/'
folder_1 = 'main_run_bcg_v2_sample/'
# folder_2 = 'less_than_10_members/'
# folder = 'image_data/'
# folder_1 = 'confirmed/'
folder_2 = 'invalid_sigma/'
already_fit = extract_fits_table(url+folder+folder_1+'v_disp_recalc_fullfitted_valid_only.fits')
#confirmed = already_fit.field('CLUS_REF')
confirmed = already_fit.field('REFNAME_1')
for i in range(0,len(confirmed)):
	conf = intstr(confirmed[i])
	#conf = confirmed[i]
	jpeg = url+folder+folder_1+conf+'.jpg'
	r_img = url+folder+folder_1+conf+'_r.fits'
	i_img = url+folder+folder_1+conf+'_i.fits'
	g_img = url+folder+folder_1+conf+'_g.fits'

	# if os.path.isfile(r_img) == True:
	# 	os.rename(r_img,url+folder+folder_1+folder_2+conf+'_r.fits')
	# if os.path.isfile(i_img) == True:
	# 	os.rename(i_img,url+folder+folder_1+folder_2+conf+'_i.fits')
	# if os.path.isfile(g_img) == True:
	# 	os.rename(g_img,url+folder+folder_1+folder_2+conf+'_g.fits')
	# if os.path.isfile(jpeg) == True:
	# 	os.rename(jpeg,url+folder+folder_1+folder_2+conf+'.jpg')
	if os.path.isdir(url+folder+conf) == True:
		move(url+folder+conf,url+folder+folder_1+conf)
################################################################################
#!$%"^^^" Nanomaggies!
url = '/home/kate/Desktop/New_BCGs/image_data/confirmed/'
folder = 'SIGMA/'
files = extract_fits_table(url+'confirmed_final_unique.fits')
clus_id = files.field('CLUS_ID')
bands = ['g','r','i']
for i in range(0,len(clus_id)):
	for g in range(0,len(bands)):
		data = fits.open(url+folder+clus_id[i]+'_'+bands[g]+'.fits')
		image = data[0].data
		header = data[0].header
		nmgy = header['NMGY']
		conversion = image/nmgy
		#hdu = fits.PrimaryHDU(conversion)
		#hdulist = fits.HDUList([hdu])
		if os.path.isfile(url+folder+clus_id[i]+'_'+bands[g]+'_adu.fits') == False:
			fits.writeto(url+folder+clus_id[i]+'_'+bands[g]+'_adu.fits',conversion,header)
		#print nmgy
		# time.sleep(5)
################################################################################
#edit inspection code for lower-ranked bcgs
# lower_ranked_candidates_for_inpection = extract_fits_table(url+'bcg_final_v2_10_z.fits')
# RA_BCG = lower_ranked_candidates_for_inpection.field('RA_BCG')
# DEC_BCG = lower_ranked_candidates_for_inpection.field('DEC_BCG')
# RA_X = lower_ranked_candidates_for_inpection.field('RA_X')
# DEC_X = lower_ranked_candidates_for_inpection.field('DEC_X')
# clus_id = lower_ranked_candidates_for_inpection.field('REFNAME_1')
# Z_ALL = lower_ranked_candidates_for_inpection.field('Z_ALL')
# RA_ALL = lower_ranked_candidates_for_inpection.field('ALLRA')
# DEC_ALL = lower_ranked_candidates_for_inpection.field('ALLDEC')
# Z_ALL = lower_ranked_candidates_for_inpection.field('ALLZ_NOQSO')
# Z_ALL_ERR = lower_ranked_candidates_for_inpection.field('ALLZ_ERR_NOQSO')
# ALLBCGFLAG = lower_ranked_candidates_for_inpection.field('ALLBCGFLAG')
# NHASZ = lower_ranked_candidates_for_inpection.field('NHASZ')
# RA_CODEX = lower_ranked_candidates_for_inpection.field('RA_2')
# DEC_CODEX = lower_ranked_candidates_for_inpection.field('DEC_2')
# clus_id = lower_ranked_candidates_for_inpection.field('CLUS_ID')
# header = "CLUS_ID,BCG_FLAG,RA_BCG,DEC_BCG,BCG_X,BCG_Y,Z_BCG,Z_ERR_BCG,CORRECT_BCG"
#bcg_data = []
# mistakes = [ '2_8885']
# for i in range(166,len(clus_id)):
# 	g = 0
# 	#bcg_data = []
# 	# candidate_bcgflag = ALLBCGFLAG[i]
# 	# sort_bcgflag = sorted(candidate_bcgflag[0:NHASZ[i]])
# 	# median_z = np.median(np.array(Z_ALL[i][0:NHASZ[i]]))
# 	codex_xy = sdss_celestial_to_x_y(url+folder+folder_1+intstr(clus_id[i])+'_r.fits', ra=RA_X[i], dec=DEC_X[i])
# 	# bcgmask = [candidate_bcgflag == sort_bcgflag[0]]
# 	# z = Z_ALL[i][bcgmask]
# 	# z_err = Z_ALL_ERR[i][bcgmask]
# 	# ra = RA_ALL[i][bcgmask]
# 	# dec = DEC_ALL[i][bcgmask]
# 	ra = RA_BCG[i]
# 	dec = DEC_BCG[i]
# 	z = Z_ALL[i]
# 	coords = sdss_celestial_to_x_y(url+folder+folder_1+intstr(clus_id[i])+'_r.fits', ra=ra, dec=dec)
# 	print 'Cluster ID = ', clus_id[i]
# 	print 'Object redshift = ', z
# 	image = mpimage.imread(url+folder+folder_1+intstr(clus_id[i])+'.jpg')
# 	plt.imshow(image)
# 	plt.plot(coords[0],1489-coords[1], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='b',markeredgewidth=1.2)
# 	plt.plot(codex_xy[0],1489-codex_xy[1], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='r',markeredgewidth=1.2)
# 	plt.axis('tight')
# 	plt.show()
# 	user_input = raw_input("Correct Centroid/Only Good Centroid? ")
# 	if user_input == "y":
# 		bcg_data = [clus_id[i],sort_bcgflag[g],ra,dec,coords[0],coords[1],z,z_err,user_input]
# 	elif user_input != "y":
# 		while user_input != "y":
# 			g = g + 1
# 			bcgmask = [candidate_bcgflag == sort_bcgflag[g]]
# 			z = Z_ALL[i][bcgmask]
# 			z_err = Z_ALL_ERR[i][bcgmask]
# 			ra = RA_ALL[i][bcgmask]
# 			dec = DEC_ALL[i][bcgmask] #check the centroid catalogues next
# 			coords = sdss_celestial_to_x_y(url+folder+clus_id[i]+'_r.fits', ra=ra, dec=dec)
# 			plt.imshow(image)
# 			plt.plot(coords[0],1489-coords[1], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='b',markeredgewidth=1.2)
# 			plt.plot(codex_xy[0],1489-codex_xy[1], marker='s',markersize=8,markerfacecolor='none',markeredgecolor='r',markeredgewidth=1.2)
# 			plt.axis('tight')
# 			plt.show()
# 			n = (len(sort_bcgflag) - (1+g))
# 			print 'Cluster ID = ', clus_id[i]
# 			print 'Number of candidates remaining: ', n
# 			print 'Object redshift = ', z
# 			print 'Median redshift =', median_z
# 			user_input = raw_input("Correct BCG? ")
# 		user_input2 = raw_input("Correct BCG? ")
# 		bcg_data = [clus_id[i],sort_bcgflag[g],ra,dec,coords[0],coords[1],z,z_err,user_input2]
# 	output = np.str(bcg_data[0])+","
# 	print output
# csvsave(url+folder+'lower_ranking_search.ascii',header,bcg_data)
################################################################################
#edit v.disp code to store all galaxies used in v disp and save them

































################################################################################

#BITS AND BOBS

################################################################################
#download jpegs/extract header data for initial checking. May need to tinker with them to get certain things to work
# url = '/home/kate/Desktop/'
# folder = 'New_BCGs/'
# new_fields = extract_fits_table(url+folder+'image_data/alt_xray/alt_xcoords_check.fits') #255
# #fields to check through
# filename = new_fields.field('filename')
# refname = new_fields.field('CLUS_ID')
# #only need inspect r-band images, if any at all. will need to match up fields to the filenames
# for i in range(0,len(refname)):
# 	image = url+folder+'image_data/alt_xray/'+filename[i]
# 	#print image
# 	camcol = getheadparam(image,param='CAMCOL',extension=0)
# 	field = getheadparam(image,param='FRAME',extension=0)
# 	rerun = getheadparam(image,param='RERUN',extension=0)
# 	run = getheadparam(image,param='RUN',extension=0)
# 	print run, rerun, camcol, field
# 	get_sdss_jpeg(rerun=rerun,run=run,camcol=camcol,field=field,location=url+folder+'image_data/alt_xray/'+refname[i]+'.jpeg')
################################################################################
#post matching. You can do this through python, however may just be easier
# filenames = extract_fits_table(file_location=url+folder+'image_data/new_sample_for_inspection.fits')
# file_ids = filenames.field('filename')
# refnames = filenames.field('CLUS_ID')
# run_array = filenames.field('RUN')
# camcol_array = filenames.field('CAMCOL')
# field_array = filenames.field('FIELD')
# for i in range(0,len(file_ids)):
# 	print file_ids[i]
# 	os.rename(url+folder+'image_data/'+file_ids[i], url+folder+'image_data/'+refnames[i]+'_r.fits')
# 	run = intstr(run_array[i])
# 	camcol = intstr(camcol_array[i])
# 	field = intstr(field_array[i])
# 	if np.int(field[i]) < 10:
# 		os.rename(url+folder+'image_data/g/frame-g-00'+run+'-'+camcol+'-000'+field+'.fits',refnames[i]+'_g.fits')
# 		os.rename(url+folder+'image_data/i/frame-i-00'+run+'-'+camcol+'-000'+field+'.fits',refnames[i]+'_i.fits')
# 	elif np.int(field[i]) >10 and np.int(field[i]) <100:
# 		os.rename(url+folder+'image_data/g/frame-g-00'+run+'-'+camcol+'-00'+field+'.fits',refnames[i]+'_g.fits')
# 		os.rename(url+folder+'image_data/i/frame-i-00'+run+'-'+camcol+'-00'+field+'.fits',refnames[i]+'_i.fits')
# 	else:
# 		os.rename(url+folder+'image_data/g/frame-g-00'+run+'-'+camcol+'-0'+field+'.fits',refnames[i]+'_g.fits')
# 		os.rename(url+folder+'image_data/i/frame-i-00'+run+'-'+camcol+'-0'+field+'.fits',refnames[i]+'_i.fits')
################################################################################
