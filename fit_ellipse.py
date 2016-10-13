
def fit_ellipse(inurl,outurl,image,mask,x0,y0,sma0,pa0,ellip0,minsma,maxsma):
	import os
	import os.path
	from pyraf import iraf
	import numpy as np
	iraf.stsdas()
	iraf.analysis()
	iraf.isophote()
	iraf.unlearn('ellipse')
	for i in range(0,3): #for three iterations
		# Define the name of the input and output file
		inputImg = inurl+image
		outBin = inputImg.replace(".fits", "_ellipse_1.bin")
		outTab = inputImg.replace(".fits", "_ellipse_1.tab")
		outCdf = inputImg.replace(".fits", "_ellipse_1.cdf")
		mask = inurl+mask
		# TODO: Check the .pl mask file, which should be
		# inputMsk = inputImg.replace(".fit", ".pl")

		# Call the STSDAS.ANALYSIS.ISOPHOTE package
		iraf.stsdas()
		iraf.analysis()
		iraf.isophote()
		iraf.unlearn('ellipse')
		# iraf.unlearn('controlpar')
		# iraf.unlearn('geompar')
		# Define parameters for the ellipse run
		# 1. Initial guess of the central X, Y (need to be as accurate as possible)
		iraf.ellipse.geompar.x0 = x0
		iraf.ellipse.geompar.y0 = y0
		# 2. Initial guess of the ellipticity and PA of the first ISOPHOTE
		#    Do not need to be very accurate, unless you want to fix them for all
		#    isophotes and only derive surface brightness
		iraf.ellipse.geompar.ellip0 = ellip0
		iraf.ellipse.geompar.pa0 = pa0
		# 3. Initial radius for ellipse fitting (The major axis length of the first
		#    elliptical isophote); Can not be too small, and can not be too large
		iraf.ellipse.geompar.sma0 = sma0
		# 4. The minimum and maximum radius for the ellipse fitting
		iraf.ellipse.geompar.minsma = minsma*sma0
		iraf.ellipse.geompar.maxsma = maxsma*sma0
		# 5. Parameters about the stepsize during the fitting.
		#    Unless you know what you what, normally should use log-stepsize instead of
		#    linear one; and step=0.05 will generate more isophotes than step=0.1, but
		#    may not help if you want a robust surface brightness profile.
		iraf.ellipse.geompar.linear = "no"
		iraf.ellipse.geompar.step = 0.1
		# 6. Do you want to allow the ellipse to decide the galaxy center during the
		#    fitting.  In general, it's a good idea to turn this on.  If the center you
		#    provide is accurate enough, ELlipse results will not deviate from it.
		iraf.ellipse.geompar.recenter = "yes"
		# 7. The next three parameters control the behavior of the fit
		#    hcenter = yes/no : Do all the isophotes have the same central X, Y?
		#    hellip  = yes/no : Do all the isophotes have the same ellipticity?
		#    hpa     = yes/no : Do all the isophotes have the same position angle?
		# Based on our experience, the formal Ellipse fitting should be done in three
		# separate runs
		#    1) hcenter=no, hellip=no, hpa=no : Give Ellipse the total freedom to fit
		#       the isophotes; And take the median/mean central X,Y from inner N
		#       isophotes, then use these X,Y as the center of the galaxy
		#    2) hcenter=yes, hellip=no, hpa=yes : Hold the central X, Y to the
		#       previously determined ones; Let the ellipticity and position angle to be
		#       free, then extract an appropriate average ellipticity and PA from this
		#       run
		#    3) hcenter=yes, hellip=yes, hpa=yes : Hold the center, and hold the
		#       ellipticity and PA to the average values decided from the previous run.
		#       Just extracted an robust surface brightness profile using the average
		#       geometry
		#print sma0, i
		# 8. Parameters about the iterations
		#    minit/maxit: minimun and maximum number of the iterations
		iraf.ellipse.controlpar.minit = 10
		iraf.ellipse.controlpar.maxit = 100
		# 9. Threshold for the object locator algorithm
		#    By lowering this value, the locator become less strict.
		iraf.ellipse.controlpar.olthresh = 1.00000
		# 10. Make sure the Interactive Mode is turned off
		iraf.ellipse.interactive="no"
		iraf.ellipse.dqf=mask
		iraf.ellipse.verbose='yes' #turms output on
		if i == 0:
			iraf.ellipse.controlpar.hcenter = "yes"
			iraf.ellipse.controlpar.hellip = "no"
			iraf.ellipse.controlpar.hpa = "yes"
			if os.path.exists(outTab+'.txt'):
				os.remove(outTab+'.txt')
			if os.path.exists(outTab):
				os.remove(outTab)
			if os.path.exists(outCdf):
				os.remove(outCdf)
			iraf.ellipse(input=inputImg, output=outTab, x0=x0,y0=y0,pa0=pa0,ellip0=ellip0)
			iraf.tables.ttools.tdump(table=outTab, datafile=outTab+'.txt', cdfile=outCdf)
			os.remove(outCdf)
		elif i == 1:
			infile = open(outTab+'.txt')
			pa = []
			x = []
			y = []
			for line in infile.readlines():
				line=line.strip()
				line=line.split()
				line[line == 'INDEF'] == np.nan
				valpa = np.float(line[7])
				valx  = np.float(line[9])
				valy = np.float(line[11])
				pa.append(valpa)
				x.append(valx)
				y.append(valy)
			pa1 = np.median(pa)
			x1 = np.median(x)
			y1 = np.median(y)
			#print pa1,x1,y1
			# iraf.ellipse.geompar.pa0 = pa1
			# iraf.ellipse.geompar.x0 = x1
			# iraf.ellipse.geompar.y0 = y1
			iraf.ellipse.controlpar.hcenter = "yes"
			iraf.ellipse.controlpar.hellip = "yes"
			iraf.ellipse.controlpar.hpa = "yes"
			if os.path.exists(outTab+'.txt'):
				os.remove(outTab+'.txt')
			if os.path.exists(outTab):
				os.remove(outTab)
			if os.path.exists(outCdf):
				os.remove(outCdf)
			iraf.ellipse(input=inputImg, output=outTab, x0=x0,y0=y0,pa0=pa1,ellip0=ellip0)
			iraf.tables.ttools.tdump(table=outTab, datafile=outTab+'.txt', cdfile=outCdf)
			os.remove(outCdf)
		elif i == 2:
			infile = open(outTab+'.txt')
			iraf.ellipse.controlpar.hcenter = "no"
			iraf.ellipse.controlpar.hellip = "no"
			iraf.ellipse.controlpar.hpa = "no"
			pa = []
			x = []
			y = []
			ellip = []
			for line in infile.readlines():
				line=line.strip()
				line=line.split()
				line[line == 'INDEF'] == np.nan
				valpa = np.float(line[7])
				valx  = np.float(line[9])
				valy = np.float(line[11])
				valellip = np.float(line[5])
				pa.append(valpa)
				x.append(valx)
				y.append(valy)
				ellip.append(valellip)
			# pa2 = np.median(pa)
			# x2 = np.median(x)
			# y2 = np.median(y)
			ellip2 = np.mean(ellip)
			if ellip2< 0.05:
				ellip2=0.05
			# iraf.ellipse.geompar.pa0 = pa1
			# iraf.ellipse.geompar.x0 = x1
			# iraf.ellipse.geompar.y0 = y1
			iraf.ellipse.geompar.ellip0 = ellip2
			if os.path.exists(outTab+'.txt'):
				os.remove(outTab+'.txt')
			if os.path.exists(outTab):
				os.remove(outTab)
			if os.path.exists(outCdf):
				os.remove(outCdf)
			iraf.ellipse(input=inputImg, output=outTab, x0=x0,y0=y0,pa0=pa1,ellip0=ellip2)
			iraf.tables.ttools.tdump(table=outTab, datafile=outTab+'.txt', cdfile=outCdf)
			os.remove(outCdf)
		# Check and remove outputs from the previous Ellipse run, or Ellipse will report
		# error (Quite stupid!)
		# if os.path.exists(outTab+'.txt'):
		# 	os.remove(outTab+'.txt')
		# if os.path.exists(outTab):
		# 	os.remove(outTab)
		# if os.path.exists(outCdf):
		# 	os.remove(outCdf)
		# Start the fitting
		# TODO: Demonstrate the direct photometry mode using input catalog
		# inBin = input_bin_file
		## The inBin is a Binary result from previous Ellipse run, and the isophote
		## stored in it will overwrite all the above settings. Ellipse will simply
		## extract surface brightness profile using these isophote instead of doing any
		## fitting
		# iraf.ellipse(input=inputImg, output=outBin, inellip=inBin)
		#print iraf.stsdas.analysis.isophote.ellipse.getParam("hcenter")
		# The Ellipse output is a binary table file, which is very hard to deal with
		#  "Dump" it into a nice ASCII table
		# iraf.tables.ttools.tdump(table=outTab, datafile=outTab+'.txt', cdfile=outCdf)
		# os.remove(outCdf)
	tablefile = outTab+'.txt'
	newTab = tablefile.replace(inurl, outurl)
	infile = outTab.replace(inurl, outurl)
	os.rename(tablefile,newTab)
	os.rename(outTab,infile)
	return ;
