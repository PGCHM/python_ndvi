#!/usr/bin/env python

# Dr. M. Disney, Aug 2011
# C. Peng, Sep 2013

import sys
import argparse
import numpy as np
from pyhdf.SD import SD, SDC
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
# local files
import hdf_utils, plot_me, linmod_me


def do_linmod():
    "do_linmod(): function to read MODIS refl HDF file and invert spectral end-members against refl data, per pixel"

    ifile = 'MOD09GA.A2004154.h19v10.005.2008152141836.hdf'
    ofile = ifile + str('.fwd')
    rmsefig = ifile + str('.rmse.png')
    paramfig = ifile + str('.params')

    # reflectance value scalar for modis data
    scale = 10000
    
    # Default multi-linear model in this example is:
    # roh(lambda) = p1 * soil(lambda) + p2 * dry_grass(lambda) + p3 * green_leaf(lambda) + p4 * offset(lambda)
    # where offset is a scaling term (effectively a sort of overall brightness). We can add (or remove) end-members
    # as we see fit (see below).

    # spectral end-member file defaults list - file format is 2 columns: lambda, roh(lambda)
    # kfiles = ['ksoil.dat.modis', 'kleaf.dat.modis']
    kfiles = ['ksoil.dat.modis', 'kdrygrass.dat.modis', 'kleaf.dat.modis', 'koffset.dat.modis']

    # default no. of spectral endmembers (params)
    nparams = len(kfiles)

    # default no. of bands - modis bands
    nbands = 7

    # default sds base name
    sdsName = 'sur_refl_b0'
    
    # default extents of image
    xmin = 0
    xmax = 250
    ymin = 0
    ymax = 250

    xlabel = 'Col'
    ylabel = 'Row'

    # specify array from nbands + nparams
    kdat = np.zeros((nparams,nbands))

    # get wb from first end-member file
    wb = np.genfromtxt(kfiles[0], unpack=True)[0]

    # read in end-members from kfiles
    for i in np.arange(nparams):
        kdat[i] = np.genfromtxt(kfiles[i], unpack=True)[1]


    # We can now plot the end-member files interactively if we want to, as we have a wb array (wb) and
    # each of the endmember files - you don't need the labels and legend, but these can be useful.....
    #for i in np.arange(nparams):
    #    plt.plot(wb, kdat[i], label='%s'%(kfiles[i]))
    #plt.xlabel('Wavelength ($\lambda$)')
    #plt.ylabel(r'Reflectance, $\mathrm{\rho}$')
    #plt.legend()
    #plt.show()
    #plt.savefig('endmembers.png')
        

    # default plot value
    plot = 1

    if options.ifile:
        ifile = options.ifile
    if options.ofile:
        ofile = options.ofile
    if options.rmsefile:
        rmsefile = options.rmsefile
    if options.rmsefig:
        rmsefig = options.rmsefig
    if options.paramfig:
        paramfig = options.paramfig
    if options.sdsName:
        sdsName = options.sdsName

    # read HDF file
    hdfip = hdf_utils.r(ifile)
    
    if options.v:
        # o/p file datasets
        sys.stderr.write('ifile: %s\n'%(ifile))
        hdfip.datasets()

    # xdim sds_red.dimensions('XDim').values()[0][0]
    # ydim sds_red.dimensions('YDim').values()[0][0]
 
    # do we have extract i.e. xmin, xmax, ymin, ymax? 
    if options.xmin:
        xmin = options.xmin    
    if options.ymin:
        ymin = options.ymin
    if options.xmax:
        xmax = options.xmax
    if options.ymax:
        ymax = options.ymax

    if options.xlabel:
        xlabel = options.xlabel

    if options.ylabel:
        ylabel = options.ylabel
 
    xrange = xmax - xmin
    yrange = ymax - ymin

    # specify input and output data arrays
    data = np.zeros((nbands,xrange,yrange)) 
    dataop = np.zeros((nbands,xrange,yrange)) 
    rmseop = np.zeros((xrange,yrange)) 
    paramsop = np.zeros((nparams,xrange,yrange)) 

    for i in np.arange(1,nbands+1):
        # sort out name of SDS
        sds = hdfip.select(sdsName + str(i) + '_1')
        # read required extract of SDS into data array AND divide by scalar
        data[i-1] = np.float32(sds[xmin:xmax,ymin:ymax])/scale
  
    # now do inversion for each pixel to calculate parameters
    for x in np.arange(xrange):
        for y in np.arange(yrange):
            # pass obs, kerns to linmod BUT remember to change order of refl data to wb order
            
            rho = modis_swap(data[:,x,y])

            # this is the line that does the inversion
	    pp = np.linalg.lstsq(np.transpose(kdat), rho)
	    
	    # returns: array containing:
	    # pp[0] = params
	    # pp[1] = sum of residuals, so rmse = np.sqrt(pp[1]/nbands)
	    # pp[2] = rank of matrix (nparams)
            # Singular values of param vector (p[0]) from inversion procedure. A cut-off for these values rcond, can be specified
            # np.linalg.lstsq(np.transpose(kdat), rho, rcond=0.01) so that values < rcond will be set to 0.
 
	    # so get params, rmse and fwd model
 	    params = pp[0]
	    rmse = np.sqrt(pp[1]/nbands)
	    fwd = np.dot(params, kdat)

            dataop[:,x,y] = fwd
            rmseop[x,y] = rmse * 100.
            paramsop[:,x,y] = params

            # keep us informed of progress per row
        if options.v:
            sys.stderr.write('\b\b\b\b\b\b(%d%%)'%(int(100.*x/xrange)))
            

    # open/create o/p files and write data to hdf datasets
    hdf_utils.w(ofile, dataop, 'fwd_model', nbands, xrange, yrange, 'reflectance', 'row', 'col')
    hdf_utils.w(rmsefig, rmseop, 'RMSE', 1, xrange, yrange, '%', 'row', 'col')
    hdf_utils.w(paramfig, paramsop, 'params', nparams, xrange, yrange, 'val.', 'row', 'col')

    # now o/p the data to plot (png files if required) - scale each for plotting i.e. mean +/- 1 sd for rmse set lo to 0 if below.
    hi, lo = hilo(rmseop, 1.)
    lo = 0 if lo < 0 else lo
    plot_me.plot(rmseop, rmsefig, xlabel, ylabel, 'RMSE %', rmsefig, lo, hi)

    hi, lo = hilo(dataop, 2.)
    lo = 0 if lo < 0 else lo
    for i in np.arange(nbands):
        plot_me.plot(dataop[i], ofile + 'b: ' + str(i), xlabel, ylabel, 'refl.', ofile + '.b' + str(i) + '.png', lo, hi)

    for i in np.arange(nparams):
        hi, lo = hilo(paramsop[i], 2)
        plot_me.plot(paramsop[i], 'param: ' + str(i), xlabel, ylabel, 'param val.', paramfig + '.' + str(i) + '.png', lo, hi)
            

def modis_swap(rho):
    '''
    modis_swap(rho)L: swap elements of reflectance array assuming modis wb order i.e. swap:
    660, 840, 485, 570, 1240, 1650, 2220
    to
    485, 570, 660, 840, 1240, 1650, 2220
    '''
    s = np.array([rho[2], rho[3], rho[0], rho[1], rho[4], rho[5], rho[6]])
    return(s)

def hilo(array, n):
    '''
    hilo(array, n): return array_mean + n*stdev, array_mean - n*stdev
    '''
    lo = array.mean() - n*np.std(array)
    hi = array.mean() + n*np.std(array)
    return(hi,lo)

def main ():
    do_linmod()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ifile", help="read in data from FILE", metavar="FILE")
    parser.add_argument("-o", "--ofile",  help="output fwd model file", metavar="FILE")
    parser.add_argument("-r", "--rmsefile",  help="rmse file", metavar="FILE")
    parser.add_argument("-p", "--paramfig", help="param fig base name: emf, eps, pdf, png, ps, raw, rgba, svg, svgz", metavar="FILE")
    parser.add_argument("-m", "--rmsefig",  help="rmse fig: emf, eps, pdf, png, ps, raw, rgba, svg, svgz", metavar="FILE")
    parser.add_argument("-s", "--sdsName",  help="SDS name")
    parser.add_argument("--xlabel", help="xlabel name")
    parser.add_argument("--ylabel", help="ylabel name")
    parser.add_argument("--xmin", type=int, help="xmin")
    parser.add_argument("--xmax", type=int, help="xmax")
    parser.add_argument("--ymin", type=int, help="ymin")
    parser.add_argument("--ymax", type=int, help="ymax")
    parser.add_argument("-v", action="store_true", help="switch verbose on")
    options = parser.parse_args()
    main()
