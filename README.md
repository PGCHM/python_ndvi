# python_ndvi
calculate and visualize NDVI indices

Execution script:

cd <dir>
  
./hdf_read_example.py -o x.png -s sur_refl_b02_1

./do_ndvi.py --xmin 1500 --xmax 1750 --ymax 1000 --ylabel "Lat" --xlabel "Lon" -t "Small NDVI" -o ndvi_small.png

./plot_pixel.py -e -x 1000 -y 1000

./do_linmod.py -v --xmin 1000 --xmax 1500 --ymin 1000 --ymax 1500
