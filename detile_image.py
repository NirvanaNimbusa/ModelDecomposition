#!/usr/bin/env python

import Image
from numpy import array,fromstring,uint8,ones,zeros,float32,arange,transpose
from numpy.fft import fft,ifft,fft2,ifft2	#rfft,irfft,rfft2,irfft2
import numpy as np

from csv import reader

#from re import compile,MULTILINE
#
#color_exp=compile(r'^0\s+!COLOU?R\s+(\w+)\s+CODE\s+(\d+)\s+VALUE\s+#([A-Fa-f0-9]{6})\s+EDGE\s+\#(\d+|[A-Fa-f0-9]{6}).*?$',MULTILINE)
#
#def parse_color(fname):
#	colors = []
#	with open(fname,'r') as fhandle:
#		fdata=fhandle.read()
#		colors=color_exp.findall(fdata)
#	return colors

epsilon = np.finfo(np.float32).eps

simple_ufunc = lambda f:np.frompyfunc(f,1,1)

@simple_ufunc
def itofcolor(color):
	return color/255.

@simple_ufunc
def ftoicolor(color):
	return int(color*255)

#@np.vectorize
#def getnth(v,n):
#	return v[n]

def nargmax(myarray, n):
	#print n,myarray
	
	ind = np.argsort(myarray)
	return ind[-n:][::-1]#, myarray[ind[-n:]]

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


#def find_nearest(difs)

def optimize_palette(palette):
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
	#print "\tInitializing floating point colors"
	#npimg=itofcolor(array(image))
	print "\tColor mapping"
	return np.apply_along_axis(find_nearest,-1,itofcolor(array(image)))

def get_vector_fetch(palette):
	#@simple_ufunc
	def fetch_color(index):
		return palette[index[0],0]
	return fetch_color
	

def depalettize(indexed,palette):
	width,height=indexed.shape
	fetch_color = get_vector_fetch(palette)
	#ret = zeros((width,height,3))
	return ftoicolor(np.apply_along_axis(fetch_color,-1,indexed[:,:,np.newaxis])).astype(uint8)
	#ret[:,:,0] = getnth(recolored[:,:],0)
	#ret[:,:,1] = getnth(recolored[:,:],1)
	#ret[:,:,2] = getnth(recolored[:,:],2)
	#print "\tCasting to int"
	##ret = ftoicolor(ret).astype(uint8)
	#ret = ftoicolor(recolored).astype(uint8)
	

def display_grids(*surfs):
	if 'Axes3D' not in globals() or 'plt' not in globals():
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
	for name,surf,(xo,xs,yo,ys) in surfs:
		fig = plt.figure() # too old to do subplots correctly!
		ax = fig.add_subplot(111, projection='3d')
		
		height,width=surf.shape
		
		x = arange(xo,width,xs)
		y= arange(yo,height,ys)
		X, Y = np.meshgrid(x, y)
		
		#print X.shape,Y.shape,surf.shape
		ax.plot_surface(X, Y, surf[yo::ys,xo::xs],cmap='spectral')
		#ax.set_zlim(*(np.mean(surf) + array([-3,3])*np.std(surf)))
		
		ax.set_xlabel('X Label')
		ax.set_ylabel('Y Label')
		ax.set_zlabel(name)

		plt.show()
	

def add_suffix(basename,suffix):
	return ".".join([
						"-".join([".".join(basename[:-1]),suffix]),
						basename[-1]])


def main(*args):
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
	print "Done"
	
	print
	print
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
	#print max_x_freq,max_y_freq
	
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
	
	#print xmaxima
	#print ymaxima
	#print
	
	#print np.logical_not(xmaxima[:,np.newaxis] % kxrange[1:]).astype(uint8)
	#print np.logical_not(ymaxima[:,np.newaxis] % kyrange[1:]).astype(uint8)
	print "Suggested tiling schemes"
	print "X:",1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(xmaxima[:,np.newaxis] % kxrange[1:])),4)
	print "Y:",1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(ymaxima[:,np.newaxis] % kyrange[1:])),4)
	#
	#
	#print 2+np.argmax(np.apply_along_axis(np.sum,0,np.logical_not(xmaxima[:,np.newaxis] % kxrange[1:])))
	#print 2+np.argmax(np.apply_along_axis(np.sum,0,np.logical_not(ymaxima[:,np.newaxis] % kyrange[1:])))
	
	#print np.apply_along_axis(np.mean,-1,xstrider(kxrange[1:])[:,np.newaxis])
	#print kxrange[::kxrange]
	#print kyrange[::kyrange]
	
	print "Done"
	
	
if __name__=="__main__":
	import sys
	main(*sys.argv[1:])
	
