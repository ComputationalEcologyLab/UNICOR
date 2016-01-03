# -------------------------------------------------------------------------------------------------------------
# DebugTools.py
# Author(s): Brian Hand
#
# originally created  : May 2011
# GOAL: Keep general debugging tools together in a utilities module
# Keywords: General Debugging, Tools, Utilities
# ------------------------------------------------------------------------------------------------------
import pdb

try:
	import sys, os
except ImportError as eMsg:
	print("ImportError: (%s) OS and SYS required!"%(eMsg))
	sys.exit(-1)
try:
	import numpy as np
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)
	
#For debugging purposes uncomment
import pylab as pylab
import matplotlib.pyplot as plt	
import matplotlib.cm as cm
	
# ------------------------------------------------------------------------------
def displayPaths(ncols,nrows,xrasloc,yrasloc,pathadd,rasterMapFile):
	'''
	displayPaths() creates a figure with the added paths over-
	layed on the resistance grids
	'''
	cmapRas = cm.gray_r
	cmapPaths = cm.hot
	
	# Turn -9999 values back to nodata 0 values
	pathadd[pathadd==-9999] = 0.0
	
	fig = plt.figure()
	rasterMap = np.loadtxt(rasterMapFile,skiprows = 6)
	rasterMap[rasterMap <0.0] = 0.0
	pathadd[pathadd <0.0] = 0.0	
	
	# For plotting, add 1 to pathaddcopy
	pathadd = pathadd + 1
	
	# Turn 1.0 values back to nodata 0 values
	pathadd[pathadd==1.0] = 0.0
	
	rasterMap_max = np.max(rasterMap)
	pathadd_max = np.max(pathadd)

	if rasterMap_max > 1:
		rasterMap = rasterMap/rasterMap_max
	rasterMap[pathadd >= 1.] = 0.0
	rasterMap = rasterMap + pathadd
	rasterMap = np.flipud(rasterMap)
	breakPoints = np.arange(0.0, 1.0, 0.1)

	if pathadd_max == 1:
		pathadd_max = pathadd_max +1
		
	# To get original path add values, subtract 1 to pathadd
	pathadd = pathadd - 1
	
	# Turn 1.0 values back to nodata 0 values
	pathadd[pathadd==-1.0] = 0.0

	breakPoints2 = np.arange(1.0,pathadd_max, 1)
	breakPoints2ticks = np.arange(1.0,pathadd_max,2)

	x = np.arange(0,ncols,1)
	y = np.arange(0,nrows,1)	
	
	cs = plt.contourf(x,y,rasterMap,breakPoints,cmap=cmapRas,  \
					alpha=1.0, clip=False)
	cs2 = plt.contourf(x,y,rasterMap,breakPoints2,cmap=cmapPaths,  \
					alpha=1.0, clip=False)
	cb = pylab.colorbar(cs,cmap=cmapRas,ticks=breakPoints,orientation='horizontal')
	cb.set_label('Resistance')
	
	cb2 = pylab.colorbar(cs2,cmap=cmapPaths,ticks=breakPoints2ticks,orientation='vertical')
	
	# Label for color bar
	cb2.set_label('Path Frequency')
	
	# Title of figure
	plt.title('Resistance Map with Paths')
	
	plt.show()	
	
	# End::displayPaths()