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
def getheadparam(file_id,extension,param):
	data = fits.open(file_id)
	variable = data[extension].header[param]
	return variable;

def intstr(input):
	output = np.str(np.int(input))
	return output;
#%%
#Extract fits and apply to images, de_vaucouleurs to start

#URL = '/home/kate/Desktop/'
#catalogue = fits.open(URL+'bcg_final_v2.fits')
#catalogue = catalogue[1].data
#ids = catalogue.field('REFNAME_1')
#xbcg_frame = catalogue.field('X_BCG_1')
#ybcg_frame = catalogue.field('Y_BCG_1')
#xbcg_frame = xbcg_frame[0]
#ybcg_frame = ybcg_frame[0]
#frame_x = 2048
#frame_y = 1489
#
#test_fit = fits.open(URL+'SIGMA/12/sigfit/12_sdss_r_modV.sigma.fits')
#test_img = fits.open(URL+'SIGMA/12_r.fits')
#NMGY = test_img[0].header['NMGY']
#IMG = np.array(test_img[0].data)/NMGY #apply nanomaggy conversion
#dimensions_y = test_fit[1].header['NAXIS2']
#dimensions_x = test_fit[1].header['NAXIS1']
#BCG_fitcoords_x = test_fit[1].header['1_XC']
#BCG_fitcoords_x = BCG_fitcoords_x.split()
#BCG_fitcoords_y = test_fit[1].header['1_YC']
#BCG_fitcoords_y = BCG_fitcoords_y.split()
#
##Defining regions of fitting box relative to entire frame
#box_edge_x_1 = np.round(xbcg_frame - np.float(BCG_fitcoords_x[0]))
#box_edge_y_1 =  np.round(ybcg_frame - np.float(BCG_fitcoords_y[0]))
#box_edge_x_2 = np.round(xbcg_frame + (dimensions_x - np.float(BCG_fitcoords_x[0])))
#box_edge_y_2 =  np.round(ybcg_frame + (dimensions_y - np.float(BCG_fitcoords_y[0])))
#
#sub_img = np.zeros([1489,2048])
#xlim_fitbox_1 = 0
#ylim_fitbox_1 = 0
#xlim_fitbox_2 = dimensions_x
#ylim_fitbox_2 = dimensions_y
#if box_edge_x_1 < 0:
#    xlim_fitbox_1= np.abs(box_edge_x_1)
#    box_edge_x_1 = 0
#elif box_edge_y_1 < 0:
#    ylim_fitbox_1 = np.abs(box_edge_y_1)
#    box_edge_y_1 = 0
#elif box_edge_x_2 > 2048:
#    xlim_fitbox_2 = dimensions_x - (box_edge_x_2 - 2048)
#    box_edge_x_2 = 2048
#elif box_edge_y_2 > 1489:
#    ylim_fitbox_2 = dimensions_x - (box_edge_y_2 - 1489)
#    box_edge_y_2 = 1489
#print box_edge_x_1,box_edge_x_2,box_edge_y_1,box_edge_y_2, xbcg_frame,ybcg_frame
#fitregion = test_fit[1].data
#sub_img[box_edge_y_1:box_edge_y_2,box_edge_x_1:box_edge_x_2] = fitregion[ylim_fitbox_1:ylim_fitbox_2,xlim_fitbox_1:xlim_fitbox_2]
##plt.imshow(IMG-sub_img,cmap='gray',norm=LogNorm()) #test works

#%%
#for real....

URL = '/home/kate/Desktop/'
catalogue = fits.open(URL+'bcg_final_v2.fits')
catalogue = catalogue[1].data
ids = catalogue.field('REFNAME_1')
ids = ids[ids!=117]#these must be the ones not fit by sigma
ids = ids[ids!=391]
ids = ids[ids!=64]
ids = ids[ids!=221]
ids = ids[ids!=238]
ids = ids[ids!=428]
ids = ids[ids!=330]
ids = ids[ids!=20]
xbcg = catalogue.field('X_BCG_1')
ybcg = catalogue.field('Y_BCG_1')
frame_x = 2048
frame_y = 1489
model = 'V'
band = 'r'
#for i in range(0,len(ids)):
#    print ids[i]
#    xbcg_frame = xbcg[i]
#    ybcg_frame = ybcg[i]
#
#    test_fit = fits.open(URL+'SIGMA/'+intstr(ids[i])+'/sigfit/'+intstr(ids[i])+'_sdss_'+band+'_mod'+model+'.sigma.fits')
#    test_img = fits.open(URL+'SIGMA/'+intstr(ids[i])+'_'+band+'.fits')
#    NMGY = test_img[0].header['NMGY']
#    IMG = np.array(test_img[0].data)/NMGY #apply nanomaggy conversion
#    dimensions_y = test_fit[1].header['NAXIS2']
#    dimensions_x = test_fit[1].header['NAXIS1']
#    BCG_fitcoords_x = test_fit[1].header['1_XC']
#    BCG_fitcoords_x = BCG_fitcoords_x.split()
#    BCG_fitcoords_y = test_fit[1].header['1_YC']
#    BCG_fitcoords_y = BCG_fitcoords_y.split()
#
#    #Defining regions of fitting box relative to entire frame
#    box_edge_x_1 = np.round(xbcg_frame - np.float(BCG_fitcoords_x[0]))
#    box_edge_y_1 =  np.round(ybcg_frame - np.float(BCG_fitcoords_y[0]))
#    box_edge_x_2 = np.round(xbcg_frame + (dimensions_x - np.float(BCG_fitcoords_x[0])))
#    box_edge_y_2 =  np.round(ybcg_frame + (dimensions_y - np.float(BCG_fitcoords_y[0])))
#
#    sub_img = np.zeros([1489,2048])
#    xlim_fitbox_1 = 0
#    ylim_fitbox_1 = 0
#    xlim_fitbox_2 = dimensions_x
#    ylim_fitbox_2 = dimensions_y
#    if box_edge_x_1 < 0:
#        xlim_fitbox_1= np.abs(box_edge_x_1)
#        box_edge_x_1 = 0
#    elif box_edge_y_1 < 0:
#        ylim_fitbox_1 = np.abs(box_edge_y_1)
#        box_edge_y_1 = 0
#    elif box_edge_x_2 > 2048:
#        xlim_fitbox_2 = dimensions_x - (box_edge_x_2 - 2048)
#        box_edge_x_2 = 2048
#    elif box_edge_y_2 > 1489:
#        ylim_fitbox_2 = dimensions_y - (box_edge_y_2 - 1489)
#        box_edge_y_2 = 1489
#        #print box_edge_x_1,box_edge_x_2,box_edge_y_1,box_edge_y_2, xbcg_frame,ybcg_frame
#    fitregion = test_fit[1].data
#    sub_img[box_edge_y_1:box_edge_y_2,box_edge_x_1:box_edge_x_2] = fitregion[ylim_fitbox_1:ylim_fitbox_2,xlim_fitbox_1:xlim_fitbox_2]
#    data_resid = IMG - sub_img
#    prihdrr = test_img[0].header
#    prihdrr['EXPTIME'] = '1' #for galfit as in flux units, you also set gain as 1
#    prihdrr['NOTE'] = 'De Vaucouleurs bestfit residual, SIGMA'
#    prihdrr.set('GAIN', '1')
#    #prihdrr.remove('BZERO')
#    #prihdrr.remove('BSCALE')
#    prihdrr.set('BUNIT', 'ADU')

    # conversiong = hdulistg[0].header['NMGY']
    # conversionr = hdulistr[0].header['NMGY']
    # conversioni = hdulisti[0].header['NMGY']
    #
    # imgg = np.array(datag)/np.float(conversiong)
    # imgr = np.array(datar)/np.float(conversionr)
    # imgi = np.array(datai)/np.float(conversioni)

    #UNCONVERTED, IN nanomaggies
    #fits.writeto(URL+'SIGMA/nmgy_residuals_'+model+'/'+intstr(ids[i])+'_'+band+'_counts_'+model+'_resid.fits', data_resid, prihdrr, clobber=True)
    #%%
#check scale size
ids = catalogue.field('REFNAME_1')
sig_flag = catalogue.field('SIG_FLAG')
scale_type = catalogue.field('SCALE_TYPE')
xbcg = catalogue.field('X_BCG_1')
ybcg = catalogue.field('Y_BCG_1')
r200_x = catalogue.field('R200c/deg')
proper_distance_arcmin = catalogue.field('PROPER_DISTANCE_PER_ARCMIN')
r200_gapper= catalogue.field('SIGMA_GAPPER')
r200_biweight = catalogue.field('SIGMA_BIWEIGHT')
scale = 0.396
valid = []
compress = []
for i in range(0,len(ids)):
	if scale_type[i] == 'PHOTOMETRIC':
		r200_pixel = (r200_x[i]*3600)/scale
		print r200_pixel
		x_lim_1 = xbcg[i]
		x_lim_2 = 2048 - xbcg[i]
		y_lim_1 = ybcg[i]
		y_lim_2 = 1489 - ybcg[i]
		if x_lim_1 < r200_pixel and x_lim_2 < r200_pixel and y_lim_1 < r200_pixel and \
			y_lim_2 < r200_pixel:
			valid.append([ids[i],'True'])
			compress.append('True')
		else:
			valid.append([ids[i],'False'])
	elif scale_type[i] == 'SPECTROSCOPIC':
		if sig_flag[i] == 'GAPPER':
			r200_pixel = ((r200_gapper[i]/proper_distance_arcmin[i])*60)/scale
			print r200_pixel
			x_lim_1 = xbcg[i]
			x_lim_2 = 2048 - xbcg[i]
			y_lim_1 = ybcg[i]
			y_lim_2 = 1489 - ybcg[i]
			if x_lim_1 < r200_pixel and x_lim_2 < r200_pixel and y_lim_1 < r200_pixel and \
			y_lim_2 < r200_pixel:
				valid.append([ids[i],'True'])
				compress.append('True')
			else:
				valid.append([ids[i],'False'])
		elif sig_flag[i] == 'BIWEIGHT':
			r200_pixel = ((r200_biweight[i]/proper_distance_arcmin[i])*60)/scale
			print r200_pixel
			x_lim_1 = xbcg[i]
			x_lim_2 = 2048 - xbcg[i]
			y_lim_1 = ybcg[i]
			y_lim_2 = 1489 - ybcg[i]
			if x_lim_1 < r200_pixel and x_lim_2 < r200_pixel and y_lim_1 < r200_pixel and \
				y_lim_2 < r200_pixel:
				valid.append([ids[i],'True'])
				compress.append('True')
			else:
				valid.append([ids[i],'False'])

print len(compress), len(valid)
################################################################################
