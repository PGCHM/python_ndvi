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

def hdf_read_example():
    "hdf_read_example(): function to read a HDF file, take an extract and plot it in a graph. Default is ifile and sdsName"

    ifile = 'C:\\Users\\Janie\\Downloads\\MOD09GA.A2004154.h19v10.005.2008152141836.hdf'
    sdsName = 'sur_refl_b01_1'
    ofile = ''

    xmin = 0
    xmax = 2000
    ymin = 0
    ymax = 2000

    # default plot value
    plot = 1


    if options.ifile:
        ifile = options.ifile
        
    if options.ofile:
        ofile = options.ofile
        
    if options.xmin:
        xmin = options.xmin
        
    if options.xmax:
         xmax = options.xmax
        
    if options.ymin:
        ymin = options.ymin
        
    if options.ymax:
         ymax = options.ymax

    if options.sdsName:
        sdsName = options.sdsName

    
    # read sds
    md = SD(ifile, SDC.READ)
    sds = md.select(sdsName)
    
    if options.v:
        # o/p file datasets
        sys.stderr.write('ifile: %s\n'%(ifile))
        md.datasets()

    if options.plot:
        # o/p file datasets
        plot = 1 

 
    if plot:
        ex = sds[xmin:xmax,ymin:ymax]
        np.clip(ex,0.,10000,out=ex)
        imshow(ex)
        plt.colorbar(drawedges="True")
        plt.title('%s'%(sdsName))
        plt.ylabel("Row")
        plt.xlabel("Col")
        if ofile:
            # will save as emf, eps, pdf, png, ps, raw, rgba, svg, svgz
            plt.savefig(ofile)
        else:
            plt.show()     



def main ():
    hdf_read_example()


# parser - this but checks to see what arguments are passed to the script and then handles them accordingly
# The resulting values are then stored in an object called 'options', which we can then access in our program
# above as options.ifile, options.ofile, ptions.xmin etc.

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ifile", help="read in data from FILE", metavar="FILE")
    parser.add_argument("-o", "--ofile", help="output file: emf, eps, pdf, png, ps, raw, rgba, svg, svgz", metavar="FILE")
    parser.add_argument("-s", "--sdsName",  help="SDS name")
    parser.add_argument("--xmin", type=int, help="xmin")
    parser.add_argument("--xmax", type=int, help="xmax")
    parser.add_argument("--ymin", type=int, help="ymin")
    parser.add_argument("--ymax", type=int, help="ymax")
    parser.add_argument("-p", "--plot", action="store_true", help="switch plotting on")
    parser.add_argument("-v", action="store_true", help="switch verbose on")
    options = parser.parse_args()
    main()



