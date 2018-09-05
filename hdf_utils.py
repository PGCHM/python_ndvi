from pyhdf.SD import SD, SDC
import numpy as np

# http://pysclint.sourceforge.net/pyhdf/example.html
# http://hdfeos.org/software/pyhdf.php
# Dr. M. Disney, Aug 2011
# C. Peng, Sep 2013

def r(ifile):
    "hdf_read.r(ifile): function to read a HDF file"

    # read hd file
    md = SD(ifile, SDC.READ)
    return(md)

def w(ofile, data, sds_base_name, nbands, xrange, yrange, units, rowdimname, coldimname):
    "hdf_utils.w(ofile, 'sds_base_name', nbands, xdim, ydim, 'units', 'rowdimname', 'coldimname')"
    "Function to create hdf file ofile and write data to it"
    # Dr. M. Disney, Sep 2011

    # create and write hdf file via SD
    dataopf = SD(ofile, SDC.WRITE | SDC.CREATE | SDC.TRUNC)
    for i in np.arange(0,nbands):  
        sds = dataopf.create(sds_base_name + str(i), SDC.FLOAT32, (xrange, yrange))
        sds.name = '' + str(i)
        sds.units = units
        sds.setfillvalue(0)
        dim1 = sds.dim(0)
        dim1.setname(rowdimname)
        dim2 = sds.dim(1)
        dim2.setname(coldimname)
        if nbands > 1:
            sds[:] = np.float32(data[i])
        else:
            sds[:] = np.float32(data)
        
    sds.endaccess()
    dataopf.end()
