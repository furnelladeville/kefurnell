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
def extract_fits_table(file_location):
	datafile = fits.open(file_location)
	datafile = datafile[1].data
	return datafile ;
def intstr(input):
	output = np.str(np.int(input))
	return output;
################################################################################
url = raw_input("File location: ")
url = '/home/kate/'+url+'/'
to_be_moved_to = raw_input("Location of files to be moved to: ")
to_be_moved_to = '/home/kate/'+to_be_moved_to+'/'
file_containing_names = raw_input("Fits file containining file ids to be moved: ")
file_containing_names = '/home/kate/'+file_containing_names
fits_table = extract_fits_table(file_containing_names)
file_ids = fits_table.field('REFNAME')
success = 0
for i in range(0,len(file_ids)):
	if os.path.isfile(url+file_ids[i]+'jpeg') == True and \
		os.path.isfile(url+file_ids[i]+'_r_adu.fits') == True and \
		os.path.isfile(url+file_ids[i]+'_g_adu.fits') == True and \
		os.path.isfile(url+file_ids[i]+'_i_adu.fits') == True:
		os.rename(url+file_ids[i]+'_r_adu.fits',to_be_moved_to+'_r_adu.fits')
		os.rename(url+file_ids[i]+'_g_adu.fits',to_be_moved_to+'_g_adu.fits')
		os.rename(url+file_ids[i]+'_i_adu.fits',to_be_moved_to+'_i_adu.fits')
		os.rename(url+file_ids[i]+'jpeg',to_be_moved_to+'.jpeg')
		success = success + 1
	else:
		print 'ERROR: One or more files missing under filename: '+file_ids[i]
print  'Files transferred: '+intstr(success)
