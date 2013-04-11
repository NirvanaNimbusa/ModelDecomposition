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

def process_edges(edges,smoothed,margin=10):
	ys,xs=edges.shape
	mm=margin-1
	mp=margin+1
	print xs,ys
	indices = transpose(np.mgrid[margin:ys-margin,margin:xs-margin],[1,2,0])
	window = np.dstack((indices+0.5,
						smoothed[mp:ys-mm,margin:xs-margin] - smoothed[mm:ys-mp,margin:xs-margin],
						smoothed[margin:ys-margin,mp:xs-mm] - smoothed[margin:ys-margin,mm:xs-mp]))
	elmct = (xs-2*margin)*(ys-2*margin)
	flatwindow = np.reshape(window,(elmct,4))
	ekeep = np.reshape(edges[margin:ys-margin,margin:xs-margin],elmct)>0
	keepers = flatwindow[ekeep]
	print "reduced edges",elmct,"to",len(keepers)
	keepers[:,2:] = transpose(np.array([np.sum(keepers[:,2:]**2,-1)**0.5,np.arctan2(keepers[:,2],keepers[:,3])+0.5*np.pi]))
	return keepers

def cluster(img,skp,sd,max_original_feature=10000):
	print "Clustering features"
	n = len(skp)
	if not n:
		return None
	kmasko = arange(n)
	kcount = ones(n)
	
	n = max(n,max_original_feature)
	
	print dir(skp[0])
	print sd.shape
	
#	int i, j, k, maski, maskj;
#	float rot_diff, size_diff, distx, disty; 
#	const double d_threshold = _param.sift_clustering_threshold;
#	const double pi = 3.1415926535897932384626433832795;
#	const double o_threshold = pi / 18.0;
#	const double dist_threshold = 48.0;
#	const double s_max = 1.5, s_min = 1.0 / s_max;
#	///////////////////////////////////////////
#	int num = _image->getFeatureNum(); if(num == 0) return;
#	unsigned char (*des)[128] = (unsigned char (*)[128])_image->getDescriptorData().data();
#	float (*loc)[5] = (float(*)[5]) _image->getLocationData().data();
#	_kmasko.resize(0);	_kmasko.resize(num, 0);
#	_kcount.resize(0);  _kcount.resize(num, 1);
#	//first step: feature clustering
#	for(i = 0; i < num; i++) _kmasko[i] = i;
#
#	// for speed purpose
#	// the original feature migh thave 100000 features
#	num = min(num, _param.max_original_feature);
#	for(i = 1 ; i < num; i++)
#	{
#		for(j = 0; j < i; j++)
#		{
#			if(_kmasko[i] == _kmasko[j]) continue;
#			////////
#			rot_diff = (float)(fmod(loc[i][4] - loc[j][4] + pi * 3.0, 2.0 * pi) - pi);
#			if(fabs(rot_diff) > o_threshold) continue;
#			///////
#			size_diff = loc[i][3] / loc[j][3];
#			if(size_diff < s_min || size_diff > s_max) continue;
#			///
#			distx = fabs(loc[i][0] -loc[j][0]);
#			disty = fabs(loc[i][1] - loc[j][1]);
#
#
#			if(max(distx, disty) > dist_threshold * max(loc[i][3], loc[j][3])) continue;
#
#			////////
#			if(MatrixUtil::dotproduct_d(des[i], des[j]) <= d_threshold)continue;
#
#
#			maski = _kmasko[i];
#			maskj = _kmasko[j];
#
#			for(k = 0; k <i; k++)
#			{
#				if(_kmasko[k] == maski) _kmasko[k] = maskj;
#			}
#			_kmasko[i] = maskj;
#			_kcount[maskj] += _kcount[maski];
#			_kcount[maski] = 0;
#		}
#	}
#	//////
#	_image_modified = 0;
#}
	
def analyze_rep(img):
	# Assume image is rectified!
	if 'cv2' not in globals() or 'cv' not in globals():
		import cv, cv2
	detector = cv2.FeatureDetector_create("SIFT")
	descriptor = cv2.DescriptorExtractor_create("SIFT")
	skp=detector.detect(img)
	skp,sd=descriptor.compute(img,skp)
	if not len(skp):
		return None
	smoothed = cv2.GaussianBlur(img,(0,0),3*max(1024,*img.shape)/2048.)
	edges = cv2.Canny(smoothed,0.3,0.8)
	pedges = process_edges(edges,smoothed)
	cluster(img,skp,sd)
	
	#Image.fromarray(img).show()
	#Image.fromarray(smoothed).show()
	#Image.fromarray(edges).show()
	#

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

	#print "Fourier Analysis"
	#four = fft2(indexed)
	##print "Displaying"
	##display_grids(("Indexed Colors",indexed,(0,1,0,1)),("Real odd",np.real(four),(1,2,1,2)),("Real even",np.real(four),(2,2,2,2)),
	##("Imag odd",np.imag(four),(1,2,1,2)),("Imag even",np.imag(four),(2,2,2,2)))
	#
	#afour = np.abs(four)**0.5
	#
	#min_x_tile = int(args[4])
	#min_y_tile = int(args[5])
	#
	#print "Min tiling resolution",min_x_tile,"x",min_y_tile
	#
	#max_y_freq, max_x_freq = indexed.shape
	#
	#xmeans = np.mean(afour,axis=0)[:1+max_x_freq/2]
	#ymeans = np.mean(afour,axis=1)[:1+max_y_freq/2]
	#
	#kxrange = arange((1+max_x_freq/2)/min_x_tile)
	#kyrange = arange((1+max_y_freq/2)/min_y_tile)
	#
	##print xmeans[kxrange]
	##print ymeans[kyrange]
	#
	#xstrider = strider(xmeans[kxrange])
	#xmmeans = xstrider(kxrange)
	#ystrider = strider(ymeans[kyrange])
	#ymmeans = ystrider(kyrange)
	#
	#xmaxima = ((np.diff(np.sign(np.diff(xmmeans))) < 0).nonzero()[0] + 1)
	#ymaxima = ((np.diff(np.sign(np.diff(ymmeans))) < 0).nonzero()[0] + 1)
	#
	#print "Suggested tiling schemes"
	#x_suggest = 1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(xmaxima[:,np.newaxis] % kxrange[1:])),2)
	#y_suggest = 1+nargmax(np.apply_along_axis(np.sum,0,np.logical_not(ymaxima[:,np.newaxis] % kyrange[1:])),2)
	#print "X:",x_suggest,xmmeans[x_suggest]
	#print "Y:",y_suggest,ymmeans[y_suggest]
	#
	#print "Best-guess scheme"
	#print "X:",x_suggest[np.argmax(xmmeans[x_suggest])]
	#print "Y:",y_suggest[np.argmax(ymmeans[y_suggest])]
	
	print "Analyzing Repetition (OpenCV)"
	analyze_rep(indexed.astype(uint8))
	
	print "Done"
	print
	
if __name__=="__main__":
	import sys
	main(*sys.argv[1:])
	

