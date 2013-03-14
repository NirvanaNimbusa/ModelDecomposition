#!/usr/bin/env python

import Image
from numpy import array,fromstring,uint8,ones,zeros,float32,arange,transpose
from numpy.fft import fft,ifft,fft2,ifft2	#rfft,irfft,rfft2,irfft2
import numpy as np

from csv import reader

epsilon = np.finfo(np.float32).eps

simple_ufunc = lambda f:np.frompyfunc(f,1,1)

@simple_ufunc
def itofcolor(color):
	return color/255.

@simple_ufunc
def ftoicolor(color):
	return int(color*255)

def nargmax(npary, n):
	return np.argsort(npary)[-n:][::-1]

def strider(npary):
	@simple_ufunc
	def stride(n):
		return np.mean(npary[::n]) if n else npary[0]
	return stride

def parse_color(color_str):
	return fromstring(color_str[1:].decode('hex'),dtype=uint8)/255.

def make_palette(fname):
	colors = -ones((512,3,3),dtype=float32)
	with open(fname,'rU') as fhandle:
		color_reader = reader(fhandle,dialect="excel")
		color_reader.next()
		for row in color_reader:
			colors[int(row[1])] = array((parse_color(row[2]),parse_color(row[3]),array((float(row[4])/255,float(row[5]),0))))
	return colors

def optimize_palette(palette):
	# Remove transparent/luminant colors,
	# because we don't have any in our inputs, for now
	global epsilon
	to_filt=transpose(array((arange(512),palette[:,0,0],palette[:,0,1],palette[:,0,2],palette[:,2,0],palette[:,2,1])))
	filt_obj=(to_filt[:,1]>=0) * (np.abs(to_filt[:,4]-1)<(100*epsilon)) * (np.abs(to_filt[:,5])<(100*epsilon))
	good_colors = to_filt[filt_obj]
	def find_nearest(color):
		difs=np.sum(np.abs(good_colors[:,1:4]-color)**2,axis=-1)**(0.5)
		return int(good_colors[np.argmin(difs)][0]+100*epsilon)
	return find_nearest

def palettize(image,palette):
	print "\tOptimizing palette (palettize)"
	find_nearest=optimize_palette(palette)
	print "\tColor mapping"
	return np.apply_along_axis(find_nearest,-1,itofcolor(array(image)))

def get_vector_fetch(palette):
	# index is actually an array of length 3
	# to make sure our output is an array n x m x 3 containing scalars
	# instead of an array of n x m containing length-3 arrays of scalars
	def fetch_color(index):
		return palette[index[0],0]
	return fetch_color
	

def depalettize(indexed,palette):
	width,height=indexed.shape
	fetch_color = get_vector_fetch(palette)
	return ftoicolor(np.apply_along_axis(fetch_color,-1,indexed[:,:,np.newaxis])).astype(uint8)

def display_grids(*surfs):
	if 'Axes3D' not in globals() or 'plt' not in globals():
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
	for name,surf,(xo,xs,yo,ys) in surfs:
		fig = plt.figure() # too old to do subplots correctly! put this outside for loop if we update
		ax = fig.add_subplot(111, projection='3d')
		
		height,width=surf.shape
		
		x = arange(xo,width,xs)
		y= arange(yo,height,ys)
		X, Y = np.meshgrid(x, y)
		
		ax.plot_surface(X, Y, surf[yo::ys,xo::xs],cmap='spectral')
		#ax.set_zlim(*(np.mean(surf) + array([-3,3])*np.std(surf)))
		
		ax.set_xlabel('X Label')
		ax.set_ylabel('Y Label')
		ax.set_zlabel(name)

		plt.show() # de-indent this if you go back to subplots
	

def add_suffix(basename,suffix):
	return ".".join([
						"-".join([".".join(basename[:-1]),suffix]),
						basename[-1]])


def main(*args):
	print "--------===",args[1],"===--------"
	print "Making palette"
	palette= make_palette(args[0])
	
	print "Optimizing palette (main)"
	find_nearest=optimize_palette(palette)
	
	print "Opening image"
	image = Image.open(args[1])
	fname = args[1]
	basename = fname.split(".")
	ndim = (int(args[2]),int(args[3]))
	
	print "Resizing image"
	nimage = image.resize(ndim,Image.ANTIALIAS)
	#nimage.show()
	nimage.save(add_suffix(basename,"pixelated"),quality=100)
	print "Palettizing image"
	indexed = palettize(array(nimage),palette)
	
	print "Depalettizing image"
	rgbalized = depalettize(indexed, palette)
	oimg=Image.fromarray(rgbalized)
	#oimg.show()
	oimg.save(add_suffix(basename,"palettized"),quality=100)

	print "Fourier Analysis"
	four = fft2(indexed)
	#print "Displaying"
	#display_grids(("Indexed Colors",indexed,(0,1,0,1)),("Real odd",np.real(four),(1,2,1,2)),("Real even",np.real(four),(2,2,2,2)),
	#("Imag odd",np.imag(four),(1,2,1,2)),("Imag even",np.imag(four),(2,2,2,2)))
	
	afour = np.abs(four)**0.5
	
	min_x_tile = int(args[4])
	min_y_tile = int(args[5])
	
	print "Min tiling resolution",min_x_tile,"x",min_y_tile
	
	max_y_freq, max_x_freq = indexed.shape
	
	xmeans = np.mean(afour,axis=0)[:1+max_x_freq/2]
	ymeans = np.mean(afour,axis=1)[:1+max_y_freq/2]
	
	kxrange = arange((1+max_x_freq/2)/min_x_tile)
	kyrange = arange((1+max_y_freq/2)/min_y_tile)
	
	#print xmeans[kxrange]
	#print ymeans[kyrange]
	
	xstrider = strider(xmeans[kxrange])
	xmmeans = xstrider(kxrange)
	ystrider = strider(ymeans[kyrange])
	ymmeans = ystrider(kyrange)
	
	xmaxima = ((np.diff(np.sign(np.diff(xmmeans))) < 0).nonzero()[0] + 1)
	ymaxima = ((np.diff(np.sign(np.diff(ymmeans))) < 0).nonzero()[0] + 1)
	
	print "Suggested tiling schemes"
	x_suggest = 1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(xmaxima[:,np.newaxis] % kxrange[1:])),2)
	y_suggest = 1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(ymaxima[:,np.newaxis] % kyrange[1:])),2)
	print "X:",x_suggest,xmmeans[x_suggest]
	print "Y:",y_suggest,ymmeans[y_suggest]
	
	print "Best-guess scheme"
	print "X:",x_suggest[np.argmax(xmmeans[x_suggest])]
	print "Y:",y_suggest[np.argmax(ymmeans[y_suggest])]
	
	
	print "Done"
	print
	
if __name__=="__main__":
	import sys
	main(*sys.argv[1:])
	
