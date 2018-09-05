# called within do_ndvi

def ndvi(red, nir):
    "ndvi(md, red, nir): module to calculate ndvi from MODIS HDF refl data md, \
    i.e. with arrays of red and nir bands passed through"
    ndvi = (nir - red) / (nir + red)
    return(ndvi)

    
