# Dr. M. Disney, Aug 2011
# C. Peng, Sep 2013

from matplotlib.pyplot import *
import matplotlib.pyplot as plt

def plot(data, title, xlabel, ylabel, colorbar_label, ofile, min, max):
	"plot(data, title, xlabel, ylabel, colorbar_label, save, ofile): module to plot data with colorbar"
	fig = plt.figure()
	ax = fig.add_subplot(111)
	cax = ax.imshow(data)
	cax.set_clim(min, max)
	ax.set_title('%s'%(title))
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	cbar = fig.colorbar(cax, orientation='vertical')
	cbar.set_label(colorbar_label)
	if ofile:
		# will save as emf, eps, pdf, png, ps, raw, rgba, svg, svgz
		plt.savefig(ofile)
	else:	
		plt.show()
	plt.close()




	
