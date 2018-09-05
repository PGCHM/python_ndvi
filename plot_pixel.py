#!/usr/bin/env python

# Dr. M. Disney, Aug 2011
# C. Peng, Sep 2013

import sys
import numpy as np
from pyhdf.SD import SD, SDC
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import matplotlib.font_manager 
from optparse import OptionParser
# local files
import hdf_utils

# def modis_swap(roh):     return(np.array([roh[2], roh[3], roh[0], roh[1], roh[4], roh[5], roh[6]]))

def plot_pixel ():
	"plot_pixel(): module to plot a pixel from two images (orig and fwd model)"
	
	wbfile = 'wb.modis.dat'
	wb = np.genfromtxt(wbfile)
	nbands = len(wb)

	scale = 10000

	kfiles = ['ksoil.dat.modis', 'kdrygrass.dat.modis', 'kleaf.dat.modis', 'koffset.dat.modis']
	nparams = len(kfiles)
	kdat = np.zeros((nparams,nbands))
	
	for i in np.arange(nparams):
		kdat[i] = np.genfromtxt(kfiles[i], unpack=True)[1]
		
	orig = 'C:\\Users\\Janie\\Downloads\\MOD09GA.A2004154.h19v10.005.2008152141836.hdf'
	# fwd = 'MOD09GA.A2004154.h19v10.005.2008152141836.hdf.fwd'
	opdat = 'op.plot.dat'
	saveplot = 'op.plot.dat.png'
	
	origbase = 'sur_refl_b0'
	
	if options.wbfile:
		wbfile = options.wbfile
		
	wb = np.genfromtxt(wbfile)
	# swap for plotting
	wb = modis_swap(wb)
	
	if options.orig:
		orig = options.orig
	if options.fwd:
		fwd = options.fwd
	if options.saveplot:
		saveplot = options.saveplot
	if options.opdat:
		opdat = options.opdat
	if options.origbase:
		origbase = options.origbase
	if options.fwdbase:
		fwdbase = options.fwdbase	

	origip = hdf_utils.r(orig)
	# fwdip = hdf_utils.r(fwd)

	# pixel location defaults
	x = 1000
	y = 300

	if options.x:
		x = options.x
	if options.y:
		y = options.y
	
	origdata = np.zeros(nbands)

	# read in pixel values
	for i in np.arange(0,nbands):
		origsds = origip.select(origbase + str(i+1) + '_1')
		origdata[i] = np.float32(origsds[x,y])/scale

	origdata = modis_swap(origdata)

	# do the inversion - np.linalg.lstsq() returns an array of 4 things, the most important of which are the first two:
	# params[0] = an array of our 4 model parameters
	# params[1] = sum of residuals i.e. the sum of the differences between original data the fwd model values
	# so  RMSE (root mean square error) = np.sqrt(params[1]/nbands)

	params = np.linalg.lstsq(np.transpose(kdat), origdata)

	# calculate our fwd model values and RMSE (multiply by 100 to get %)
	fwddata = np.dot(params[0],kdat)
	rmse = np.sqrt(params[1]/nbands) * 100.

	# set up a plot figure, and a subplot ax that we can modify
	fig = plt.figure()
	ax = plt.subplot(111)

	# plot the original data and fwd modelled values
	ax.plot(wb, origdata, 'ro', label=r'$\rho_{orig}$')
	ax.plot(wb, fwddata, 'k+', label=r'$\rho_{fwd}$')

	# if we want them, plot our end member spectra too, to compare
	if options.e:
		for i in np.arange(nparams):
			ax.plot(wb, kdat[i], label='%s'%(kfiles[i]))

	#box = ax.get_position()
	#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])


	# plotting niceties. Set the y axis limits, x and y axis labels, a title and then add some
	# text to the plot showing the RMSE and the parameter values
	ax.set_ylim(0, 0.9)
	ax.set_xlabel('$\lambda$ (nm)')
	ax.set_ylabel(r'Reflectance, $\rho$')
	ax.set_title('Simple linear mixture model')
	ax.text(1250, 0.1, 'RMSE (%%) = %f'%(rmse))
	ax.text(1250, 0.05, 'params = %.2f, %.2f, %.2f, %.2f'%(params[0][0], params[0][1], params[0][2], params[0][3]))

	# set the font size and then scale the legend accordingly, and put it at 0.5, 0.95 (in relative plot coordinates)
	# i.e. 50% along (x axis) and 5% from the top (y axis)
	prop = matplotlib.font_manager.FontProperties(size=10) 
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=3, fancybox=True, shadow=True, prop=prop)

	# do we want to save to file or plot to the screen?
	if options.saveplot:
		plt.savefig(saveplot)
	else:
		plt.show()
	
	# plt.close(1)

	if options.opdat:
		tmp = np.zeros((nbands, 3))
		for i in np.arange(0,nbands):
			tmp[i, 0] = wb[i]
			tmp[i, 1] = origdata[i]
			tmp[i, 2] = fwddata[i]
			# print "%i %f %f"%(wb[i], origdata[i], fwddata[i])
		np.savetxt(opdat, tmp)

def modis_swap(roh):
    # swap elements of reflectance array assuming modis wb order i.e. swap:
    # 660, 840, 485, 570, 1240, 1650, 2220
    # to
    # 485, 570, 660, 840, 1240, 1650, 2220
    return(np.array([roh[2], roh[3], roh[0], roh[1], roh[4], roh[5], roh[6]]))

def main ():
    plot_pixel()
    
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-o", "--original", dest="orig", help="original MODIS reflectance file", metavar="FILE")
    parser.add_option("-f", "--forward", dest="fwd", help="fwd model  reflectance file", metavar="FILE")    
    parser.add_option("-w", "--wavebands", dest="wbfile", help="File containing MODIS wavebands", metavar="FILE")
    parser.add_option("-p", "--plot", dest="saveplot", help="output plotfile: emf, eps, pdf, png, ps, raw, rgba, svg, svgz", metavar="FILE")
    parser.add_option("-d", "--data", dest="opdat", help="output orig & fwd model data", metavar="FILE")
    parser.add_option("--fwdbase", dest="fwdbase", help="fwd model SDS name")
    parser.add_option("--origbase", dest="origbase", help="orig refl SDS name")
    parser.add_option("--xmin", type="int", dest="xmin", help="xmin")
    parser.add_option("--xmax", type="int", dest="xmax", help="xmax")
    parser.add_option("--ymin", type="int", dest="ymin", help="ymin")
    parser.add_option("--ymax", type="int", dest="ymax", help="ymax")
    parser.add_option("-x", type="int", dest="x", help="x location")
    parser.add_option("-y", type="int", dest="y", help="y location")
    parser.add_option("-e", action="store_true", help="plot end members")
    (options, args) = parser.parse_args()
    main()


