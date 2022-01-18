# ---------------------------------------------------------------------------------------------------
# UNICORMaps.py
# Author(s): Erin L Landguth, Brian Hand, Joe Glassy
# Created: September 2009
# v 0.5 -- 2010 August 30: Code Cleanup and function additions to outputs.
# ---------------------------------------------------------------------------------------------------
# Messy, fix later
msgVerbose = False

# Import Modules with Except/Try statements

# For ELL debugging
import pdb

# Numpy functions
try:
	import numpy as np
except ImportError:
	raise ImportError("Numpy required.")


# Scipy functions
try:
	import scipy as sp
	import scipy.signal	
	scipyAvail = True
except ImportError:
	raise ImportError("Scipy required for plotting.")
	scipyAvail = False
	
	
# Timing functions and misc functions
try:
	import time, datetime, copy					
except ImportError:
	raise ImportError("Time and Datetime required.")	

# -------------------------------------------------------------
def triweight_kern(size, sizey=None):
	'''
	Returns a normalized 2D Triweight kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Biweight has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = (35./32)*((1.-x[xval][yval]**2)**3)
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = (35./32)*((1.-y[xval][yval]**2)**3)
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::triweight_kernal

# -------------------------------------------------------------
def cosine_kern(size, sizey=None):
	'''
	Returns a normalized 2D Cosine kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Biweight has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = (np.pi/4.)*np.cos((np.pi/2.)*x[xval][yval])
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = (np.pi/4.)*np.cos((np.pi/2.)*y[xval][yval])
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::cosine_kern

# -------------------------------------------------------------
def biweight_kern(size, sizey=None):
	'''
	Returns a normalized 2D Biweight kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Biweight has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = (15./16)*((1.-x[xval][yval]**2)**2)
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = (15./16)*((1.-y[xval][yval]**2)**2)
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::biweight_kern

# -------------------------------------------------------------
def triangle_kern(size, sizey=None):
	'''
	Returns a normalized 2D Trianlge kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Triangle has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = 1. - np.abs(x[xval][yval])
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = 1. - np.abs(y[xval][yval])
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::triangle_kern

# -------------------------------------------------------------
def uniform_kern(size, sizey=None):
	'''
	Returns a normalized 2D Uniform kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Uniform has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = (1./2)
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = (1./2)
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::uniform_kern

# -------------------------------------------------------------
def epanechnikov_kern(size, sizey=None):
	'''
	Returns a normalized 2D Epanechnikov kernel array for convolutions
	'''
	# For irregular grid buffer set size
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	
	# Make a integer buffer grid of x,y coords
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.ones((size*2+1,sizey*2+1),dtype=np.float)
	
	# Perform function: Epanechnikov has split case for if abs(x,y) above/below 1
	
	# Loop through grid row values
	for xval in range(int(x.shape[1])):
		
		# Loop through grid column values
		for yval in range(int(x.shape[0])):
			
			# For abs values <= 1 at x spot
			if np.abs(x[xval][yval]) <= 1:
			
				# x grid value
				xspot = (3./4)*((1.-x[xval][yval]**2))
				
			# For abs values > 1 at x spot
			else:
				
				# x grid value
				xspot = 0.0				
				
			# For abs values <= 1 at y spot
			if np.abs(y[xval][yval]) <= 1:
			
				# y grid value
				yspot = (3./4)*((1.-y[xval][yval]**2))
			
			# For abs values > 1 at y spot
			else:
				
				# y grid value
				yspot = 0.0	
		
			g[xval][yval] = xspot + yspot
	
	return g / g.sum(), g

	# End::epanechnikov_kern

# -------------------------------------------------------------
def gauss_kern(size, sizey=None):
	'''
	Returns a normalized 2D gauss kernel array for convolutions
	'''
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
	g = np.exp(-(x**2/float(size) + y**2/float(sizey)))
	return g / g.sum(), g

	# End::gauss_kern

if scipyAvail:
	# -------------------------------------------------------------
	def blur_image(im,n,xrasloc,yrasloc,KernelFunction,ny=None) :
		''' 
		blurs the image by convolving with a gaussian kernel of typical
		size n. The optional keyword argument ny allows for a different
		size in the y direction.
		'''
		
		# For KDE_Function choice Gaussian w/o optimal bandwidth
		if KernelFunction == 'Gaussian':
			gnorm,g = gauss_kern(n, sizey=ny)
			#print gnorm		
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# For KDE_Function choice Epanechnikov
		elif KernelFunction == 'Epanechnikov':
			gnorm,g = epanechnikov_kern(n, sizey=ny)
			#print gnorm	
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# For KDE_Function choice Uniform
		elif KernelFunction == 'Uniform':
			gnorm,g = uniform_kern(n, sizey=ny)
			#print gnorm
			improc = sp.signal.convolve(im, gnorm, mode='same')
		   
		# For KDE_Function choice Triangle
		elif KernelFunction == 'Triangle':
			gnorm,g = triangle_kern(n, sizey=ny)
			#print gnorm
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# For KDE_Function choice Biweight
		elif KernelFunction == 'Biweight':
			gnorm,g = biweight_kern(n, sizey=ny)
			#print gnorm
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# For KDE_Function choice Cosine
		elif KernelFunction == 'Cosine':
			gnorm,g = cosine_kern(n, sizey=ny)
			#print gnorm
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# For KDE_Function choice Cosine
		elif KernelFunction == 'Triweight':
			gnorm,g = triweight_kern(n, sizey=ny)
			#print gnorm
			improc = sp.signal.convolve(im, gnorm, mode='same')
			
		# Return smoothed image
		return improc 
		
		# End::blur_image()

# -----------------------------------------------------------------------------------------------------------
def getlevels(ncols,nrows,Levels,buffervalues)	:
	'''
	getlevels() returns an array with catagorical cutoff values.
	'''
	
	# Get maximum value from buffervalues
	buffermax = np.max(buffervalues)
	
	# Get catagorical cutoff percentage
	cutoffpercent = 100./Levels
	
	# Get catagorical cutoff buffered values in list
	cutofflevels = []
	
	# Loop through number of levels for cutoff mark
	for icut in range(Levels-1):
	
		# Get buffered values cutoff
		cutofflevels.append((cutoffpercent*(icut+1)/100)*buffermax)
		
	# Add 0 and max values
	cutofflevels.insert(0,0)
	cutofflevels.append(buffermax)
		
	# Then go through buffervalues re-assigning catagorical levels
	for icol in range(ncols):
		for irow in range(nrows):
			for icutlevel in range(Levels):
			
				# Then check cases
				if buffervalues[irow][icol] > cutofflevels[icutlevel] and buffervalues[irow][icol] <= cutofflevels[icutlevel+1]:
					
					# And assign a value + 1
					temp = icutlevel + 1
					
				else:
				
					# Or give a 0 value
					temp = buffervalues[irow][icol]
								
			# Then assign the buffervalue
			buffervalues[irow][icol] = temp
					
	return buffervalues
	
	# End::getlevels()

# -----------------------------------------------------------------------------------------------------------
def createLevels(ncols,nrows,xrasloc,yrasloc,rasterMapFile,buffervalues,Levels):
	'''
	displayLevels() creates catagorical levels of the buffered path.
	'''
	
	# Turn -9999 buffer values back to nodata 0 values
	buffervalues[buffervalues==-9999] = 0.0
	
	buffervalues[buffervalues <0.0] = 0.0
	
	# Deep copy buffervalues
	bufferlevels = copy.deepcopy(buffervalues)
				
	# Call function to get levels
	bufferlevels = getlevels(ncols,nrows,Levels,bufferlevels)
	
	# Turn 0.0 pathadd values back to nodata -9999 values
	bufferlevels[bufferlevels==0.0] = -9999
	
	return bufferlevels
	
	# End::displayLevels()

# -----------------------------------------------------------------------------------------------------------
def createBuffers(ncols,nrows,xrasloc,yrasloc,pathadd,rasterMapFile,KDE_GridSize,KernelFunction):
	'''
	displayBuffers() creates a figure with the KDE buffered paths over-
	layed on the resistance grid.
	'''
		
	# Turn -9999 values back to nodata 0 values
	pathadd[pathadd==-9999] = 0.0
	
	pathadd[pathadd <0.0] = 0.0
	
	# Deep copy values ???
	pathaddcopy = copy.deepcopy(pathadd)
	
	# For kernel density estimation
	pathaddcopy = blur_image(pathaddcopy,KDE_GridSize,xrasloc,yrasloc,KernelFunction)
	
	# Turn -1.0 vallues back to nodata -9999 values
	pathaddcopy[pathaddcopy==0.0] = -9999
	
	return pathaddcopy
	
	# End::displayBuffers()
