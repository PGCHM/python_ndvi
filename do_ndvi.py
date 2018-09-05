#!/usr/bin/env python

# http://pysclint.sourceforge.net/pyhdf/example.html
# Dr. M. Disney, Aug 2011
# C. Peng, Sep 2013

import sys
import argparse
import numpy as np
from pyhdf.SD import SD, SDC
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
# local files
import ndvi, hdf_utils, plot_me


def do_ndvi():
    "do_ndvi(): function to read MODIS refl HDF file and calculate NDVI"

    ifile = 'MOD09GA.A2004154.h19v10.005.2008152141836.hdf'
    title = ifile
    ofile = 0
    red = 'sur_refl_b01_1'
    nir = 'sur_refl_b02_1'
    
    xmin = 0
    xmax = 500
    ymin = 0
    ymax = 500

    xlabel = 'Col'
    ylabel = 'Row'
    

    # min and max values for colourbar on plot
    ndvimin = 0.1
    ndvimax = 0.8
    
    # default plot value
    plot = 1

    if options.ifile:
        ifile = options.ifile

    if options.ofile:
        ofile = options.ofile

    if options.sdsName:
        sdsName = options.sdsName

    # read HDF file
    md = hdf_utils.r(ifile)
    
    if options.v:
        # o/p file datasets
        sys.stderr.write('ifile: %s\n'%(ifile))
        md.datasets()

    sds_red = md.select(red)
    sds_nir = md.select(nir)


    # do we have extract i.e. xmin, xmax, ymin, ymax? 
    if options.xmin:
        xmin = options.xmin    
    if options.ymin:
        ymin = options.ymin
    if options.xmax:
        xmax = options.xmax
    if options.ymax:
        ymax = options.ymax

    # scale ndvi values?
    if options.ndvimin:
        ndvimin = options.ndvimin
    if options.ndvimax:
        ndvimax = options.ndvimax

    if options.xlabel:
        xlabel = options.xlabel

    if options.ylabel:
        ylabel = options.ylabel
        
    if options.title:
        title = options.title
        
    # check dims
    if (sds_red.dimensions('XDim').values()[0][0] != sds_nir.dimensions('XDim').values()[0][0]) or (sds_red.dimensions('YDim').values()[0][0] != sds_nir.dimensions('YDim').values()[0][0]):
        sys.stderr.write('%s: dimension error - x, y dims of SDS %s and %s do not match\n'%(ifile,red,nir))


    # get extract if required, and cast to float32
    if (options.xmin or options.ymin or options.xmax or options.ymax):
        sds_red = np.float32(sds_red[xmin:xmax,ymin:ymax])
        sds_nir = np.float32(sds_nir[xmin:xmax,ymin:ymax])
        n = np.zeros([xmax-xmin,ymax-ymin])
    else:
        xdim = sds_red.dimensions('XDim').values()[0][0]
        ydim = sds_red.dimensions('YDim').values()[0][0]
        sds_red = np.float32(sds_red[0:xdim, 0:ydim])
        sds_nir = np.float32(sds_nir[0:xdim, 0:ydim])
        n = np.zeros([xdim, ydim])

    # calculate ndvi from bands 1, 2
    n = ndvi.ndvi(sds_red, sds_nir)

    # clip if required
    np.clip(n,-1.,1.,out=n)
    
    if options.v:
        # o/p file datasets
        sys.stderr.write('ifile: %s\n'%(ifile))
        md.datasets()

    if options.plot:
        # o/p file datasets
        plot = 1 

    if plot:
        plot_me.plot(n, title, xlabel, ylabel, 'NDVI', ofile, ndvimin, ndvimax)

def main ():
    do_ndvi()
    

# parser - note how this is virtually the same as for hdf_read_example.py - reduce, reuse, recycle!

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ifile", help="read in data from FILE", metavar="FILE")
    parser.add_argument("-o", "--ofile", help="output file: emf, eps, pdf, png, ps, raw, rgba, svg, svgz", metavar="FILE")
    parser.add_argument("-s", "--sdsName",  help="SDS name")
    parser.add_argument("-t", "--title",  help="title")
    parser.add_argument("--xlabel", help="xlabel name")
    parser.add_argument("--ylabel", help="ylabel name")
    parser.add_argument("--xmin", type=int, help="xmin")
    parser.add_argument("--xmax", type=int, help="xmax")
    parser.add_argument("--ymin", type=int, help="ymin")
    parser.add_argument("--ymax", type=int, help="ymax")
    parser.add_argument("--ndvimin", type=float, help="ndvimin")
    parser.add_argument("--ndvimax", type=float, help="ndvimax")
    parser.add_argument("-p", "--plot", action="store_true", help="switch plotting on")
    parser.add_argument("-v", action="store_true", help="switch verbose on")
    options = parser.parse_args()
    main()
