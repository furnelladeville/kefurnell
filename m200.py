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
#import itertools
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

################################################################################
#Functions

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

def comoving_mpc(z1,z2): #z2 > z1
	Omega_l = 0.7
	Omega_m = 0.3
	c = 3*10**8
	H_0 = 70
	da_dr = lambda z: (c/H_0)*(Omega_m*(1 + z)**3.0 + Omega_l)**(-0.5) #comoving distance
	da1, err = quad(da_dr, 0, z1,limit=100)
	da2, err = quad(da_dr, 0, z2, limit=100)
	comoving = (da2/1000.0) - (da1/1000.0)
	return comoving;

# def sigma(v): #velocity dispersion of cluster using gapper estimate
# 	v_sorted = sorted(v)
# 	n_v = len(v_sorted)
# 	if n_v == 1:
# 		sigma = v[0]
# 	else:
# 		product = []
# 		coeff = (np.sqrt(np.pi))/(n_v*(n_v - 1.0))
# 		for i in range(0,len(v_sorted)-2):
# 			j = i+1
# 			gap = v_sorted[j] - v_sorted[i]
# 			weight = j*(n_v - j)
# 			product.append(weight*gap)
# 		sigma = coeff*(sum(product))
# 	return sigma;

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
			#print len(v_3), len(sigma_cut)
		if len(sigma_cut) < limit:
			gapper = np.nan
			biweight = np.nan
			flag = 0
			n_obj = len(sigma_cut)
		else:
			#GAPPER (CORRECT THIS TIME, WE HOPE)
			n_obj = len(sigma_cut)
			#BIWEIGHT (FOR COMPARISON)
			biweight = asts.biweight_location(sigma_cut,M=np.mean(sigma_cut))
			flag = 1
			#print biweight,  gapper
			print biweight, z_phot
	return biweight, sigma_cut,flag;

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
#%%
################################################################################
#Let's calculate some physical parameters
Omega_l = 0.7
Omega_m = 0.3
c = 3*10**8
H_0 = 70 #h100 in Finn at al., 2005 is left floating
c2 = (3.0)*10**5.0 #(for mass/gapper use c2)
# URL = ['/home/kate/Desktop/SIGMA-old/']
# URL2 = ['/home/kate/Desktop/']
# params = fits.open(URL[0]+'cluster_statistics_2016-04-15-DR13.fits')
# params = params[1].data
# file_id = params.field('CLUS_ID')
# n_z = params.field('NHASZ')
# z_noqso = params.field('ALLZ_NOQSO')
# z_clus = params.field('Z_LAMBDA')

# m_clus = []
# for g in range(0,len(n_z)):
# 	objs = z_noqso[g,0:np.int(n_z[g])] #gets spec_z that aren't mad for the cluster
# 	#v_objs = c2*objs
# 	v_objs_comoving = []
# 	for q in range(0,len(objs)):
# 		comoving = cosmo.comoving_distance(objs[q])
# 		comoving = comoving.value
# 		v_objs_comoving.append(H_0*comoving)
# 	sig_clus = gapper_biweight(v_objs_comoving,errors=True)
# 	sig_clus_err =  gapper_biweight_err(v_objs_comoving,sig_clus[1],5000)
# 	if sig_clus[1] == 1: #Taken from Rose Finn's 2005 paper - she recommends using X-ray gas mass as a proxy if cluster is not relaxed
# 		r200 = (1.73)*((sig_clus[0])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + z_clus[g])**3.0))**(-1.0))
# 		#print r200
# 		m200 = (1.0/0.7)*(1.2*10**15)*((sig_clus[0]/1000.0)**3)*((np.sqrt(Omega_l + Omega_m*(1 + z_clus[g])**3.0))**(-1.0))
# 		#print np.str(m200)
# 		m_clus.append([file_id[g],sig_clus[0],sig_clus[1],sig_clus[2],sig_clus[3],sig_clus_err[0],sig_clus_err[1],r200,m200])
# 	else:
# 		m_clus.append([file_id[g],sig_clus[0],sig_clus[1],sig_clus[2],sig_clus[3],sig_clus_err[0],sig_clus_err[1],r200,m200])
# filehead = ['CLUS_ID','SIGMA_GAPPER','SIGMA_BIWEIGHT','SIGMA_FLAG','N_VALID','KS_P_MEMZ','GAPPER_ERR','BIWEIGHT_ERR','R200_SIGMA_GAPPER','M200_SIGMA_GAPPER']
# csvsave(URL2[0]+'gapper_results_DR13.csv',filehead,m_clus)

URL = ['/home/kate/Desktop/']
params = fits.open(URL[0]+'bcg_final_v2.fits')
params = params[1].data
file_id = params.field('CLUS_ID_1a')
n_z = params.field('NHASZ')
z_noqso = params.field('ALLZ_NOQSO')
ra_n = params.field('ALLRA')
dec_n = params.field('ALLDEC')
ra_clus = params.field('RA_OPT')
dec_clus = params.field('DEC_OPT')
z_clus = params.field('Z_LAMBDA')
r200c = params.field('R200c/deg')
mag_r = params.field('GAL_sdss_r_modS_MAGTOT')
mag_g = params.field('GAL_sdss_g_modS_MAGTOT')
mag_i = params.field('GAL_sdss_i_modS_MAGTOT')
m_clus = []
m_clus_2 = []
# pvals = []
# n_mem = []
#cluster_velocities = []
#start_time = time.time()
#for g in range(0,len(n_z)):
#	objs = z_noqso[g,0:np.int(n_z[g])] #gets spec_z that aren't mad for the cluster
#	ra = ra_n[g,0:np.int(n_z[g])]
#	dec = dec_n[g,0:np.int(n_z[g])]
#	#v_objs = c2*objs
#	v_objs_comoving = []
#	z_obj = []
#	#l = 0
#	for q in range(0,len(objs)):
#		sep = np.arccos(cos(90-dec[q])*cos(90-dec_clus[g]) + sin(90-dec[q])*sin(90-dec_clus[g])*cos(ra[q]-ra_clus[g]))
#		if sep < r200c[g]:
#			comoving = cosmo.comoving_distance(objs[q])
#			comoving = comoving.value
#			v_objs_comoving.append(H_0*comoving)
#			z_obj.append(objs[q])
#	zbwt = biweight_z(v_objs_comoving,z_obj,z_clus[g])
#	zbwt_err = biweight_err_z(z_obj,10000,zbwt[2],zbwt[0])
#	v2 = ((np.array(zbwt[1]) - zbwt[0])/(1+zbwt[0]))*c2
#	sig_clus = gapper_biweight(v_objs_comoving)
#	sig_clus2 = gapper_biweight(v2)
#	#print sig_comov[1] - sig_clus2[1]
#	#print sig_clus
#	#time.sleep(5)
#	sig_clus_err =  gapper_biweight_err(v=v_objs_comoving,flag=sig_clus[2],n=10000, biweight_true=sig_clus[1],gapper_true=sig_clus[0])
#	sig_clus_err2 = gapper_biweight_err(v=v2,flag=sig_clus2[2],n=10000, biweight_true=sig_clus2[1],gapper_true=sig_clus2[0])
#	# print sig_clus_err
#	# time.sleep(5)
#	if sig_clus[2] == 1: #Taken from Rose Finn's 2005 paper - she recommends using X-ray gas mass as a proxy if cluster is not relaxed
#		r200 = (1.73)*((sig_clus[1])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
#		#print r200
#		m200 = (1.0/0.7)*(1.2*10**15.0)*((sig_clus[1]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
#		#print np.str(m200)
#		# cluster_velocities.append(sig_clus[4])
#		# pvals.append(sig_clus[3])
#		# n_mem.append(sig_clus[2])
#		r2002 = (1.73)*((sig_clus2[1])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
#		#print r200
#		m2002 = (1.0/0.7)*(1.2*10**15.0)*((sig_clus2[1]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[0])**3.0))**(-1.0))
#		#print sig_clus[1],sig_clus2[1],sig_clus[0],sig_clus2[0]
#		m_clus.append([file_id[g],sig_clus[0],sig_clus[1],sig_clus[2],sig_clus[3],sig_clus[5],sig_clus_err[0],sig_clus_err[1],sig_clus_err[2],sig_clus_err[3],sig_clus_err[4],sig_clus_err[5],r200,m200,zbwt[0],zbwt_err[0],zbwt_err[1],zbwt_err[2]])
#		m_clus_2.append([file_id[g],sig_clus2[0],sig_clus2[1],sig_clus2[2],sig_clus2[3],sig_clus2[5],sig_clus_err2[0],sig_clus_err2[1],sig_clus_err2[2],sig_clus_err2[3],sig_clus_err2[4],sig_clus_err2[5],r2002,m2002])
#	else:
#		m_clus.append([file_id[g],sig_clus[0],sig_clus[1],sig_clus[2],sig_clus[3],sig_clus[5],sig_clus_err[0],sig_clus_err[1],sig_clus_err[2],sig_clus_err[3],sig_clus_err[4],sig_clus_err[5],np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
#		m_clus_2.append([file_id[g],sig_clus2[0],sig_clus2[1],sig_clus2[2],sig_clus2[3],sig_clus2[5],sig_clus_err2[0],sig_clus_err2[1],sig_clus_err2[2],sig_clus_err2[3],sig_clus_err2[4],sig_clus_err2[5],np.nan,np.nan])
#	print g
#	print("--- %s seconds ---" % (time.time() - start_time))
#filehead = ['CLUS_ID','SIGMA_GAPPER','SIGMA_BIWEIGHT','SIGMA_FLAG','N_VALID','SIG_FLAG','BIWEIGHT_ERR_UPPER','BIWEIGHT_ERR_LOWER','GAPPER_ERR_UPPER','GAPPER_ERR_LOWER','BIWEIGHT_SIGMA','GAPPER_SIGMA','R200_SIGMA_BIWEIGHT','M200_SIGMA_BIWEIGHT','ZBWT','ZBWT_ERR_UPPER','ZBWT_ERR_LOWER','ZBWT_SIGMA']
#csvsave(URL[0]+'gapper_results_R200c_vcomoving.csv',filehead,m_clus)
#
#filehead = ['CLUS_ID','SIGMA_GAPPER','SIGMA_BIWEIGHT','SIGMA_FLAG','N_VALID','SIG_FLAG','BIWEIGHT_ERR_UPPER','BIWEIGHT_ERR_LOWER','GAPPER_ERR_UPPER','GAPPER_ERR_LOWER','BIWEIGHT_SIGMA','GAPPER_SIGMA','R200_SIGMA_BIWEIGHT','M200_SIGMA_BIWEIGHT']
#csvsave(URL[0]+'gapper_results_R200c_valt.csv',filehead,m_clus_2)
#%%
#Calculate gapper m200
bcgdata = fits.open('/home/kate/Desktop/De_Vaucouleurs_plus_sigma.fits')
bcgdata = bcgdata[1].data
refname = bcgdata.field('REFNAME_1')
sigma_gapper = bcgdata.field('SIGMA_GAPPER')
zbwt = bcgdata.field('ZBWT')
m_clus = []
for i in range(0,len(refname)):
    r200 = (1.73)*((sigma_gapper[i])/1000.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[i])**3.0))**(-1.0))
    m200 = (1.0/0.7)*(1.2*10**15.0)*((sigma_gapper[i]/1000.0)**3.0)*((np.sqrt(Omega_l + Omega_m*(1 + zbwt[i])**3.0))**(-1.0))
    m_clus.append([refname[i],m200,r200])
#filehead = ['REFNAME','M200_SIGMA_GAPPER','R200_SIGMA_GAPPER']
#csvsave(URL[0]+'gapper_m200_r200.csv',filehead,m_clus)
#%%
################################################################################
#Abs_Mag_MAG_GRI
# absmags = []
# for i in range(0,len(file_id)):
# 	lum_mpc = cosmo.luminosity_distance(z_clus[i]) #bcg z
# 	lum_mpc = lum_mpc.value
# 	comoving = cosmo.comoving_distance(z_clus[i])
# 	comoving = comoving.value
# 	scale = cosmo.kpc_proper_per_arcmin(z_clus[i])
# 	scale = scale.value
# 	magabs_g_s = mag_g[i] - 5*(np.log10(lum_mpc*(10**6))-1)
# 	magabs_r_s = mag_r[i] - 5*(np.log10(lum_mpc*(10**6))-1)
# 	magabs_i_s = mag_i[i] - 5*(np.log10(lum_mpc*10**6)-1)
# 	#print comoving, lum_mpc, scale
# 	absmags.append([file_id[i],lum_mpc,comoving,scale,magabs_g_s,magabs_r_s,magabs_i_s])

# filehead = ['CLUS_ID','DL_MPC','DC_MPC','PROPER_KPC_ARCMIN','M_S_g','M_S_r','M_S_i']
# csvsave(URL[0]+'absmags_n1_distmeasures.csv',filehead,absmags)
#%%
#RE DeVaucouleurs vs Cluster Mass
# URL = '/home/kate/Desktop/'
# bcgdata = fits.open(URL+'v_disp_sample_25516.fits')
# bcgdata = bcgdata[1].data
# refname = bcgdata.field('REFNAME_1')
# zbwt = np.array(bcgdata.field('ZBWT'))
# sigma_flag = bcgdata.field('SIG_FLAG')
# Re_kpc = np.array(bcgdata.field('RE_r_V_KPC'))
# clusmass_gapper = bcgdata.field('M200_SIGMA_GAPPER')
# clusmass_bwt = bcgdata.field('M200_SIGMA_BIWEIGHT')
#
# times = np.array([1.0,2.0,3.0,4.0,5.0])
# agenow = cosmo.age(0)
# agenow = agenow.value
# z_bins = []
# for i in range(0,len(times)+1):
#     if i == 0:
#         z_lookback = 0.0
#         z_bins.append(z_lookback)
#     else:
#         j = i - 1
#         t = agenow - times[j]
#         z_lookback = z_at_value(cosmo.age, t * u.Gyr)
#         z_bins.append(z_lookback)
# z_bins = np.array(z_bins)
# #print z_bins
# #Select best velocity dispersion as gapper or biweight
# gapper = (sigma_flag == 'GAPPER')
# m200_gapper = clusmass_gapper[gapper]
# re_gapper = Re_kpc[gapper]
# zbwt_gapper = zbwt[gapper]
# biweight = (sigma_flag == 'BIWEIGHT')
# m200_biweight = clusmass_bwt[biweight]
# re_biweight = Re_kpc[biweight]
# zbwt_biweight = zbwt[biweight]
# z = []
# z.extend(zbwt_gapper)
# z.extend(zbwt_biweight)
# z = np.array(z)
# Re = []
# Re.extend(re_gapper)
# Re.extend(re_biweight)
# Re = np.array(Re)
# Re = np.log10(Re)
# m200 = []
# m200.extend(m200_gapper)
# m200.extend(m200_biweight)
# m200 = np.array(m200)
# #xticks = [1,2,5,10,15,20,30,40]
# Re_bins = asts.histogram(Re,bins='knuth')
# #Re_bins = xticks
# width = Re_bins[0]
# Re_bins = Re_bins[1]
# m200_bins = []
# for i in range(0,len(Re_bins)-1):
#     j = i+1
#     msk = (Re>Re_bins[i])&(Re<=Re_bins[j])
#     m200binvals = m200[msk]
#     med_m200 = np.median(m200binvals)
#     m200_bins.append(med_m200)
# lastbin_msk = (Re>Re_bins[len(Re_bins)-1])
# m200_last = np.array(m200[lastbin_msk])
# m200_bins = np.array(m200_bins)
# m200_ext = []
# m200_ext.extend(m200_bins)
# m200_ext.append(np.median(m200_last))
# Re_bincentres = []
# for i in range(0,len(Re_bins)-1):
#     centre = (Re_bins[i]+Re_bins[j])/2
#     Re_bincentres.append(centre)
# lastbin = Re_bins[len(Re_bins)-1] + width/2
# Re_bincentres.append(lastbin)
# #print len(m200_ext),len(Re_bins)
# color=iter(cm.rainbow(np.linspace(0,1,len(z_bins)-1)))
# bin_strings = ['> 1','1 - 2','2 - 3','3 - 4','4 - 5']
# xticks = [2,3,5,10,20,40]
# xticksstring = (str(w) for w in xticks)
# xticks = np.log10(np.array(xticks))
# yticks = [1*10**13,2*10**13,5*10**13,10*10**13,20*10**13,50*10**13,100*10**13,200*10**13]
# yticksstring = ['1e13','2e13','5e13','1e14','2e14','5e14','1e15','3e15']
# for i in range(0,len(z_bins)-1):
#     j = i+1
#     msk = (z>z_bins[i])&(z<=z_bins[j])
#     zbinvals = z[msk]
#     m200binvals=m200[msk]
#     Rebinvals=Re[msk]
#     c=next(color)
#     plt.plot(Rebinvals,m200binvals,markeredgecolor=c,marker='o', linestyle='None', label=bin_strings[i]+'Gyr',markerfacecolor='none')
# plt.gca().set_yscale("log")
# plt.gca().set_xscale("log")
# plt.xlabel('Effective De Vaucouleurs Radius (kpc)')
# plt.ylabel('Cluster Mass (M_{solar})')
# plt.plot(Re_bins[3:len(Re_bins)],m200_ext[3:len(m200_ext)],color='black',marker='o', linestyle='solid',label='Median')
# plt.gca().yaxis.set_ticks(yticks)
# plt.gca().xaxis.set_ticks(xticks)
# plt.ylim([0,3*10**15])
# plt.xlim([0,150])
# plt.xlim([np.min(Re)+0.4,np.max(Re)])
# plt.gca().set_xticklabels(xticksstring)
# plt.gca().set_yticklabels(yticksstring)
# plt.legend(loc='upper left')
# plt.show()

#%%
################################################################################
#Quartile gaussian
# percentiles = [25.0,50.0,75.0]
# quartile_values = []
# gappervals = []
# rowrefs = []
# for i in range(0,len(m_clus)):
# 	if np.isnan(m_clus[i][1]) != True:
# 		gappervals.append(m_clus[i][1])
# 		rowrefs.append(i)
# #print gappervals
# for i in range(0,len(percentiles)):
# 	limval = np.percentile(gappervals, percentiles[i])
# 	quartile_values.append(limval)
# #print quartile_values
# low_25=[]
# low_25_rows=[]
# mid_50=[]
# mid_50_rows=[]
# high_75=[]
# high_75_rows=[]
# for i in range(0,len(gappervals)):
# 	if gappervals[i] < quartile_values[0]:
# 		low_25.append(gappervals[i])
# 		low_25_rows.append(rowrefs[i])
# 	elif gappervals[i] > quartile_values[0] and gappervals[i] < quartile_values[2]:
# 		mid_50.append(gappervals[i])
# 		mid_50_rows.append(rowrefs[i])
# 	elif gappervals[i] > quartile_values[2]:
# 		high_75.append(gappervals[i])
# 		high_75_rows.append(rowrefs[i])
#
#
# bins_25 = asts.knuth_bin_width(low_25, return_bins=True)
# bins_25 = bins_25[1]
# # binnew = []
# # for i in range(0,len(bins_25[1])-1):
# # 	g = i + 1
# # 	bincent = (bins_25[1][i] + bins_25[1][g])/2
# # 	binnew.append(bincent)
# # bins_25 = binnew
# bins_50 = asts.knuth_bin_width(mid_50, return_bins=True)
# bins_50 = bins_50[1]
# # binnew = []
# # for i in range(0,len(bins_50[1])-1):
# # 	g = i + 1
# # 	bincent = (bins_50[1][i] + bins_50[1][g])/2
# # 	binnew.append(bincent)
# # bins_50 = binnew
# bins_75 = asts.knuth_bin_width(high_75, return_bins=True)
# bins_75 = bins_75[1]
# # binnew = []
# # for i in range(0,len(bins_75[1])-1):
# # 	g = i + 1
# # 	bincent = (bins_75[1][i] + bins_75[1][g])/2
# # 	binnew.append(bincent)
# # bins_75 = binnew
#
# #print bins_75,bins_50,bins_25
# norm_25 = mlab.normpdf(sorted(low_25),np.mean(low_25),np.std(low_25))
# norm_50 = mlab.normpdf(sorted(mid_50),np.mean(mid_50),np.std(mid_50))
# norm_75 = mlab.normpdf(sorted(high_75),np.mean(high_75),np.std(high_75))
# from matplotlib.pyplot import subplot
# subplot(2,2,1)
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
#
# plt.hist(low_25, bins=bins_25,histtype='step',stacked=True,color='green',normed=True)
# plt.plot(sorted(low_25),norm_25,'g--')
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(low_25),np.max(low_25)])
# # plt.ylim([0,np.max(norm_25)])
# plt.title('Lower 25%, Overall Velocity Dispersion Distribution')
# plt.ylabel('Count')
# plt.xlabel('Velocity Dispersion')
# plt.axis('tight')
# #plt.axis('tight')
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# subplot(2,2,2)
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
#
# plt.hist(mid_50, bins=bins_50,histtype='step',stacked=True,color='blue',normed=True)
# plt.plot(sorted(mid_50),norm_50,'b--')
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(mid_50),np.max(mid_50)])
# # plt.ylim([0,np.max(norm_50)])
# plt.title('Central 50%, Overall Velocity Dispersion Distribution')
# plt.ylabel('Count')
# plt.xlabel('Velocity Dispersion')
# plt.axis('tight')
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# subplot(2,2,3)
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
#
# plt.hist(high_75, bins=bins_75,histtype='step',stacked=True,color='red',normed=True)
# plt.plot(sorted(high_75),norm_75,'r--')
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(high_75),np.max(high_75)])
# # plt.ylim([0,np.max(norm_75)])
# plt.title('Upper 75%, Overall Velocity DIspersion Distribution')
# plt.ylabel('Count')
# plt.xlabel('Velocity Dispersion')
# plt.axis('tight')
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# subplot(2,2,4)
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
# #plt.plot(low_25,norm_25,'g--')
# plt.hist(low_25, bins=bins_25,histtype='step',color='green',normed=True)
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(low_25),np.max(low_25)])
# # plt.ylim([0,np.max(norm_25)])
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
# #plt.plot(mid_50,norm_50,'b--')
# plt.hist(mid_50, bins=bins_50,histtype='step',color='blue',normed=True)
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(mid_50),np.max(mid_50)])
# # plt.ylim([0,np.max(norm_50)])
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# # xticks = [min(samp[i]),max(samp[i])]
# # xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# # yticks = [0,max(normclus)]
# # yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
# #plt.plot(high_75,norm_75,'r--')
# plt.hist(high_75, bins=bins_75,histtype='step',color='red',normed=True)
# # plt.gca().xaxis.set_ticks(xticks)
# # plt.gca().yaxis.set_ticks(yticks)
# # plt.xlim([np.min(high_75),np.max(high_75)])
# # plt.ylim([0,np.max(norm_75)])
# plt.title('Stacked Quartiles')
# plt.ylabel('Normed Count')
# plt.xlabel('Velocity Dispersion')
# #plt.axis('tight')
# # plt.gca().set_xticklabels(xticksstring)
# # plt.gca().set_yticklabels(yticksstring)
# #plt.axis('tight')
# #plt.show()
# ################################################################################
#Median binned sigma values
# URL = ['/home/kate/Desktop/']
#
# sigma_cat = fits.open(URL[0]+'bcg_final_v2.fits')
# sigma_cat_me = sigma_cat[1].data
# lx = sigma_cat_me.field('Lx/Ez')
# sigma_g = sigma_cat_me.field('SIGMA_GAPPER')
# sigma_b = sigma_cat_me.field('SIGMA_BIWEIGHT')
# #bins = [2*10**43,3*10**43,5*10**43,8*10**43,20*10**43,50*10**43]
# lx_clip = []
# for i in range(0,len(lx)):
# 	lxi = np.log10(lx[i])
# 	lx_clip.append(lxi)
# lx = lx_clip
# bins = asts.scott_bin_width(lx,return_bins=True)
# bins = bins[1]
# #print bins[1],bins[0]
# binmed_lx = []
# new_bin = []
# for i in range(0,len(bins)):
# 	g = i+1
# 	mid_save = []
# 	if g < len(bins):
# 		binmin = bins[i]
# 		binmax = bins[g]
# 		for p in range(0,len(lx)):
# 			if np.float(lx[p]) > binmin and np.float(lx[p]) < binmax and sigma_g[p] < 1500:
# 				mid_save.append(sigma_g[p])
# 		# print mid_save
# 		median = np.median(mid_save)
# 		newbin = (np.float(binmin) + np.float(binmax))/2.0
# 		#print newbin
# 		binmed_lx.append(median)
# 		new_bin.append(newbin)
#
# # #print binmed_lx
# #print len(new_bin),len(binmed_lx)
# new_bin = new_bin[3:len(new_bin)-1]
# binmed_lx = binmed_lx[3:len(binmed_lx)-1]
# #print len(new_bin),len(binmed_lx)
# plt.plot(new_bin,binmed_lx,'k.',markersize=12)
# plt.plot(new_bin,binmed_lx,'k-')
# # plt.plot(new_binn,binmed_lxn,'b.',markersize=12)
# # plt.plot(new_binn,binmed_lxn,'b-')
# plt.plot(lx,sigma_g,'rx')
# # plt.plot(lxn,sigman,'bx')
# #plt.gca().set_xscale("log")
# plt.gca().set_yscale("log")
# plt.ylabel('Velocity dispersion')
# plt.xlabel('log(Lx)')
# yticks = [200,400,600,800,1000,1500]
# yticksstring=['200','400','600','800','1000','1500']
# plt.gca().yaxis.set_ticks(yticks)
# plt.gca().set_yticklabels(yticksstring)
# plt.axis('tight')
# plt.xlim([42.9,45])
# plt.show()
################################################################################

################################################################################
#Plots gaussian (x150)
# subset = 25
# remain = 17
# count = 0
# n = 1
# samp = []
# ngals=[]
# psamp = []
# for g in range(0,len(cluster_velocities)):
# 	samp.append(cluster_velocities[count])
# 	psamp.append(pvals[count])
# 	ngals.append(n_mem[count])
# 	if count == subset*n and n < 7:
# 		n = n + 1
# 		for i in range(0,subset):
# 			q = i + 1
# 			mean = np.mean(samp[i])
# 			std = np.std(samp[i])
# 			#axis_lim = (-5*std,5*std,100)
# 			#axis_lim = np.linspace(min(samp[i]),max(samp[i]),100)
# 			plt.subplot(5,5,q)
# 			xticks = [min(samp[i]),max(samp[i])]
# 			xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# 			normclus = mlab.normpdf(samp[i],mean,std)
# 			yticks = [0,max(normclus)]
# 			yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
# 			plt.plot(samp[i],normclus,'r--')
# 			plt.hist(samp[i], bins=np.linspace(min(samp[i]), max(samp[i]), np.round(len(samp)/3)),histtype='stepfilled',stacked=True,color='green',normed=True)
# 			plt.gca().xaxis.set_ticks(xticks)
# 			plt.gca().yaxis.set_ticks(yticks)
# 			plt.xlim([min(samp[i]),max(samp[i])])
# 			plt.ylim([0,max(normclus)])
# 			plt.title('(p,n) = ('+np.str(np.round(psamp[i],3))+','+np.str(np.round(ngals[i],3))+')')
# 			plt.gca().set_xticklabels(xticksstring)
# 			plt.gca().set_yticklabels(yticksstring)
# 			#plt.axis('tight')
# 		plt.show()
# 		samp = []
# 		psamp = []
# 		ngals = []
# 	elif count == len(cluster_velocities):
# 		for i in range(0,remain): #USE asts.scott_bin_width(data,return_bins=true)
# 			q = i + 1
# 			mean = np.mean(samp[i])
# 			std = np.std(samp[i])
# 			#axis_lim = (-5*std,5*std,100)
# 			#axis_lim = np.linspace(min(samp[i]),max(samp[i]),100)
# 			plt.subplot(5,5,q)
# 			xticks = [min(samp[i]),max(samp[i])]
# 			xticksstring = [np.str(np.round(np.log10(min(samp[i])),3)),np.str(np.round(np.log10(max(samp[i])),3))]
# 			normclus = mlab.normpdf(samp[i],mean,std)
# 			yticks = [0,max(normclus)]
# 			yticksstring = ['0',np.str(np.round(np.log10(max(normclus)),3))]
# 			plt.plot(samp[i],normclus,'r--')
# 			plt.hist(samp[i], bins=np.linspace(min(samp[i]), max(samp[i]), np.round(len(samp)/3)),histtype='stepfilled',stacked=True,color='green',normed=True)
# 			plt.gca().xaxis.set_ticks(xticks)
# 			plt.gca().yaxis.set_ticks(yticks)
# 			plt.xlim([min(samp[i]),max(samp[i])])
# 			plt.ylim([0,max(normclus)])
# 			plt.title('(p,n) = ('+np.str(np.round(psamp[i],3))+','+np.str(np.round(ngals[i],3))+')')
# 			plt.gca().set_xticklabels(xticksstring)
# 		plt.show()
# 		samp = []
# 		psamp = []
# 		ngals = []
# 	count = count + 1
