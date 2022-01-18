# ---------------------------------------------------------------------------------------------------
# UNICORFun.py
# Author(s): Erin L Landguth, Brian Hand, Joe Glassy
# revised: 2022-01-18 ELL
# revised: 2010-12-24 BKH
# revised: 2010-11-16T14:31:01MST jmg
# revised: 2010-11-16T13:07:01MST jmg
# revised: September 2009
# revised: v0.5 2010-08-30T12:00:00MST 
#  made code Cleanup and function additions to outputs.
# ------------------------------------------------------------------------------------------------------
# Messy, fix later
msgVerbose = False

# Import Modules with Except/Try statements

# For ELL debugging
import pdb

# Main python functions
import math

# Numpy functions
try:
	import numpy as np					
except ImportError as eMsg:
	print(("ImportError: (%s) Numpy required."%(eMsg)))
	
# Timing functions
try:
	import time, datetime					
except ImportError as eMsg:
	print(("ImportError: (%s) Time and Datetime required."%(eMsg)))
	
# Platform and system functions
try:
	import os, sys                    
except ImportError as eMsg:
	print(("ImportError: (%s) OS and SYS required."%(eMsg)))

#File absolute paths for importing functions
UTILITIES_PATH =  "../utilities/"

#Import the package specific folders
FILEIO_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"FileIO"))

if FILEIO_folder not in sys.path:
     sys.path.insert(0, FILEIO_folder)
	 
# Utility Functions
try:
	import FileIO
except ImportError as eMsg:
	print(("ImportError: FileIO.py is required."%(eMsg))) 
	
# Dictionary library needed for Dijkstra's algorithm
try:
	from priodict import priorityDictionary	
	priodictAvail = True
except ImportError as eMsg:
	print(("ImportError: (%s) Priodict required."%(eMsg)))
	priodictAvail = False

try:
	from scipy.spatial import KDTree
	scipyAvail = True
except ImportError as eMsg:
	print(("ImportError: (%s) Scipy required."%(eMsg)))
	scipyAvail = False
	
# For parallel processing
try:
	import queue
	from multiprocessing import Process, JoinableQueue        
except ImportError as eMsg:
	print(("ImportError: (%s) Multiprocessing required."%(eMsg))) 



def info(title):
	'''
	info()
	Prints out information on function called,
	module called, parent process, child process,
	and pid.
	'''
	print(title)
	print('module name:', __name__)
	if os.name != 'nt':
		print('parent process:', os.getppid())
	print(("info: process id: %d"%(os.getpid())))
	# End::info()

def setFloat(strPoint):
	isource_x = float(strPoint.split('_')[0])
	isource_y = float(strPoint.split('_')[1])
	return (isource_x,isource_y)

	
def setStr(xxx_todo_changeme):
	(float_x,float_y) = xxx_todo_changeme
	return str(float_x) + '_' + str(float_y)	

def approxWhere(chkAry, searchValue, TOLERANCE = 10E-5):
	criterion = (chkAry >= searchValue - TOLERANCE) & \
	(chkAry <= searchValue + TOLERANCE) 
	indices = np.where(criterion)
	return indices
	
if priodictAvail:
	# --------------------------------------------------------------------------
	def Dijkstra(nWeightPairs,xvalues,yvalues,start, endpts=None):
		"""
		Find shortest paths from the start vertex to all
		vertices nearer than or equal to the end.

		The input graph G is assumed to have the following
		representation: A vertex can be any object that can
		be used as an index into a dictionary.  G is a
		dictionary, indexed by vertices.  For any vertex v,
		G[v] is itself a dictionary, indexed by the neighbors
		of v.  For any edge v->w, G[v][w] is the length of
		the edge.  This is related to the representation in
		<http://www.python.org/doc/essays/graphs.html>
		where Guido van Rossum suggests representing graphs
		as dictionaries mapping vertices to lists of neighbors,
		however dictionaries of edges have many advantages
		over lists: they can store extra information (here,
		the lengths), they support fast existence tests,
		and they allow easy modification of the graph by edge
		insertion and removal.  Such modifications are not
		needed here but are important in other graph algorithms.
		Since dictionaries obey iterator protocol, a graph
		represented as described here could be handed without
		modification to an algorithm using Guido's representation.

		Of course, G and G[v] need not be Python dict objects;
		they can be any other object that obeys dict protocol,
		for instance a wrapper in which vertices are URLs
		and a call to G[v] loads the web page and finds its links.
		
		The output is a pair (D,P) where D[v] is the distance
		from start to v and P[v] is the predecessor of v along
		the shortest path from s to v.
		
		Dijkstra's algorithm is only guaranteed to work correctly
		when all edge lengths are positive. This code does not
		verify this property for all edges (only the edges seen
		before the end vertex is reached), but will correctly
		compute shortest paths even for some graphs with negative
		edges, and will raise an exception if it discovers that
		a negative edge has caused it to make a mistake.
		"""	
		D = {}	# dictionary of final distances
		P = {}	# dictionary of predecessors
		Q = priorityDictionary()   # est.dist. of non-final vert.
		Q[start] = 0
		vOrder = []
		paths = []
		pathlens = []
		storex = []	
		storey = []
		
		for v in Q:
			D[v] = Q[v]
			(vx,vy) = setFloat(v)
			
			xIndex = approxWhere(xvalues,vx,10E-3)
			yIndex = approxWhere(yvalues,vy,10E-3)
			
			x = xIndex[0][0]
			y = yIndex[0][0]
			storex.append(x)
			storey.append(y)
			
			if (endpts.count(v) > 0):
				endpts.remove(v)
				Path = shortestPath(D,P,start,v)
				vOrder.append(v)
				paths.append(Path)
				pathlens.append(D[v])
			if  len(endpts) == 0: break

			# w is a neighboring vertex
			for w in nWeightPairs[y][x]:
				
				if str(w[0]) == 'nan':
					continue
				
				wStr = setStr((w[0],w[1]))
				vwLength = D[v] + w[2]
				
				if wStr in D:	
					if vwLength < D[wStr]:
						raise ValueError("Dijkstra: found better path to already-final vertex")
				
				elif wStr not in Q or vwLength < Q[wStr]:
					Q[wStr] = vwLength
					P[wStr] = v
		
		return (paths, vOrder,pathlens)
		
		# End::Dijkstra()
	
	# --------------------------------------------------------------------------
	def Dijkstra_cond(nWeightPairs,xvalues,yvalues,start, endpts=None):
		"""
		Find shortest paths from the start vertex to all
		vertices nearer than or equal to the end.

		The input graph G is assumed to have the following
		representation: A vertex can be any object that can
		be used as an index into a dictionary.  G is a
		dictionary, indexed by vertices.  For any vertex v,
		G[v] is itself a dictionary, indexed by the neighbors
		of v.  For any edge v->w, G[v][w] is the length of
		the edge.  This is related to the representation in
		<http://www.python.org/doc/essays/graphs.html>
		where Guido van Rossum suggests representing graphs
		as dictionaries mapping vertices to lists of neighbors,
		however dictionaries of edges have many advantages
		over lists: they can store extra information (here,
		the lengths), they support fast existence tests,
		and they allow easy modification of the graph by edge
		insertion and removal.  Such modifications are not
		needed here but are important in other graph algorithms.
		Since dictionaries obey iterator protocol, a graph
		represented as described here could be handed without
		modification to an algorithm using Guido's representation.

		Of course, G and G[v] need not be Python dict objects;
		they can be any other object that obeys dict protocol,
		for instance a wrapper in which vertices are URLs
		and a call to G[v] loads the web page and finds its links.
		
		The output is a pair (D,P) where D[v] is the distance
		from start to v and P[v] is the predecessor of v along
		the shortest path from s to v.
		
		Dijkstra's algorithm is only guaranteed to work correctly
		when all edge lengths are positive. This code does not
		verify this property for all edges (only the edges seen
		before the end vertex is reached), but will correctly
		compute shortest paths even for some graphs with negative
		edges, and will raise an exception if it discovers that
		a negative edge has caused it to make a mistake.
		"""	
		D = {}	# dictionary of final distances
		D_cond = {}	# dictionary of final distances with cond
		P = {}	# dictionary of predecessors
		Q = priorityDictionary()   # est.dist. of non-final vert.
		Q_cond = priorityDictionary()   # est.dist. of non-final vert.
		Q[start] = 0
		Q_cond[start] = 1
		vOrder = []
		paths = []
		pathlens = []
		pathcond = []
		storex = []	
		storey = []
		
		for v in Q:
		
			D[v] = Q[v]
			D_cond[v] = Q_cond[v]
			(vx,vy) = setFloat(v)
			
			xIndex = approxWhere(xvalues,vx,10E-3)
			yIndex = approxWhere(yvalues,vy,10E-3)
			
			x = xIndex[0][0]
			y = yIndex[0][0]
			storex.append(x)
			storey.append(y)
			
			if (endpts.count(v) > 0):
				endpts.remove(v)
				Path = shortestPath(D,P,start,v)
				vOrder.append(v)
				paths.append(Path)
				pathlens.append(D[v])
				pathcond.append(D_cond[v])
			if  len(endpts) == 0: break

			# w is a neighboring vertex
			for w in nWeightPairs[y][x]:
				
				if str(w[0]) == 'nan':
					continue
				
				wStr = setStr((w[0],w[1]))
				vwConduct = D_cond[v] * w[3]
				vwLength = D[v] + w[2]
				
				if wStr in D:	
					if vwLength < D[wStr]:
						pdb.set_trace()
						raise ValueError("Dijkstra: found better path to already-final vertex.")
				
				elif wStr not in Q or vwLength < Q[wStr]:
				#if wStr not in Q:
					Q[wStr] = vwLength
					P[wStr] = v
					Q_cond[wStr] = vwConduct
		
		return (paths, vOrder,pathlens,pathcond)
		
		# End::Dijkstra_cond()

	# --------------------------------------------------------------------------
	def Dijkstra_wind(nWeightPairs,xvalues,yvalues,start, endpts=None):
		"""
		Find shortest paths from the start vertex to all
		vertices nearer than or equal to the end.

		The input graph G is assumed to have the following
		representation: A vertex can be any object that can
		be used as an index into a dictionary.  G is a
		dictionary, indexed by vertices.  For any vertex v,
		G[v] is itself a dictionary, indexed by the neighbors
		of v.  For any edge v->w, G[v][w] is the length of
		the edge.  This is related to the representation in
		<http://www.python.org/doc/essays/graphs.html>
		where Guido van Rossum suggests representing graphs
		as dictionaries mapping vertices to lists of neighbors,
		however dictionaries of edges have many advantages
		over lists: they can store extra information (here,
		the lengths), they support fast existence tests,
		and they allow easy modification of the graph by edge
		insertion and removal.  Such modifications are not
		needed here but are important in other graph algorithms.
		Since dictionaries obey iterator protocol, a graph
		represented as described here could be handed without
		modification to an algorithm using Guido's representation.

		Of course, G and G[v] need not be Python dict objects;
		they can be any other object that obeys dict protocol,
		for instance a wrapper in which vertices are URLs
		and a call to G[v] loads the web page and finds its links.
		
		The output is a pair (D,P) where D[v] is the distance
		from start to v and P[v] is the predecessor of v along
		the shortest path from s to v.
		
		Dijkstra's algorithm is only guaranteed to work correctly
		when all edge lengths are positive. This code does not
		verify this property for all edges (only the edges seen
		before the end vertex is reached), but will correctly
		compute shortest paths even for some graphs with negative
		edges, and will raise an exception if it discovers that
		a negative edge has caused it to make a mistake.
		"""	
		D = {}	# dictionary of final distances
		D_wind = {}	# dictionary of final distances with cond
		P = {}	# dictionary of predecessors
		Q = priorityDictionary()   # est.dist. of non-final vert.
		Q_wind = priorityDictionary()   # est.dist. of non-final vert.
		Q[start] = 0
		Q_wind[start] = 1
		vOrder = []
		paths = []
		pathlens = []
		pathwind = []
		storex = []	
		storey = []
		
		for v in Q:
		
			D[v] = Q[v]
			D_wind[v] = Q_wind[v]
			(vx,vy) = setFloat(v)
			
			xIndex = approxWhere(xvalues,vx,10E-3)
			yIndex = approxWhere(yvalues,vy,10E-3)
			
			x = xIndex[0][0]
			y = yIndex[0][0]
			storex.append(x)
			storey.append(y)
			
			if (endpts.count(v) > 0):
				endpts.remove(v)
				Path = shortestPath(D,P,start,v)
				vOrder.append(v)
				paths.append(Path)
				pathlens.append(D[v])
				pathwind.append(D_wind[v])
			if  len(endpts) == 0: break

			# w is a neighboring vertex
			for w in nWeightPairs[y][x]:
				
				if str(w[0]) == 'nan':
					continue
				
				wStr = setStr((w[0],w[1]))
				#vwWind = D_wind[v] * (1. - w[2]) # Take 1 - to get back to original prob
				vwWind = D_wind[v] + w[4] # This is the cost
				#vwLength = D[v] + w[2]
				vwLength = D[v] + w[2] # This is the cost to search for path
				#vwCost = D_cost[v] + w[4]
				
				if wStr in D:	
					if vwLength < D[wStr]:
						raise ValueError("Dijkstra: found better path to already-final vertex.")
				
				elif wStr not in Q or vwLength < Q[wStr]:
				#if wStr not in Q:
					Q[wStr] = vwLength
					P[wStr] = v
					Q_wind[wStr] = vwWind
			
		return (paths, vOrder,pathlens,pathwind)
		
		# End::Dijkstra_wind()

	
	# --------------------------------------------------------------------------			
	def shortestPath(D,P,start,end):
		'''
		shortestPath()
		Find a single shortest path from the given start vertex
		to the given end vertex.
		The input has the same conventions as Dijkstra().
		The output is a list of the vertices in order along
		the shortest path.
		vOrder is the order in which the vertices were found
		'''
		(startX,startY) = setFloat(start)
		
		Path = []
		if end not in P:
			Path.append(None)
		
		elif end in P:
			while 1:
				Path.append(end)
				(endX,endY) = setFloat(end)
				if (isApprox(endX,startX,10E-3) and isApprox(endY,startY,10E-3)): 
					break
				end = P[end]
			Path.reverse()
		return Path
		
		# End::shortestPath()
	
	# --------------------------------------------------------------------------
	def Dijkstrathreshold(edgethreshold,nWeightPairs,xvalues,\
			yvalues,start,endpts=None):
		"""
		Find shortest paths from the start vertex to all
		vertices nearer than or equal to the end.
		
		Except if path is beyond edgethreshold length.
		Then return string 'Path beyond edge threshold'.

		The input graph G is assumed to have the following
		representation: A vertex can be any object that can
		be used as an index into a dictionary.  G is a
		dictionary, indexed by vertices.  For any vertex v,
		G[v] is itself a dictionary, indexed by the neighbors
		of v.  For any edge v->w, G[v][w] is the length of
		the edge.  This is related to the representation in
		<http://www.python.org/doc/essays/graphs.html>
		where Guido van Rossum suggests representing graphs
		as dictionaries mapping vertices to lists of neighbors,
		however dictionaries of edges have many advantages
		over lists: they can store extra information (here,
		the lengths), they support fast existence tests,
		and they allow easy modification of the graph by edge
		insertion and removal.  Such modifications are not
		needed here but are important in other graph algorithms.
		Since dictionaries obey iterator protocol, a graph
		represented as described here could be handed without
		modification to an algorithm using Guido's representation.

		Of course, G and G[v] need not be Python dict objects;
		they can be any other object that obeys dict protocol,
		for instance a wrapper in which vertices are URLs
		and a call to G[v] loads the web page and finds its links.
		
		The output is a pair (D,P) where D[v] is the distance
		from start to v and P[v] is the predecessor of v along
		the shortest path from s to v.
		
		Dijkstra's algorithm is only guaranteed to work correctly
		when all edge lengths are positive. This code does not
		verify this property for all edges (only the edges seen
		before the end vertex is reached), but will correctly
		compute shortest paths even for some graphs with negative
		edges, and will raise an exception if it discovers that
		a negative edge has caused it to make a mistake.
		"""

		D = {}	# dictionary of final distances
		P = {}	# dictionary of predecessors
		Q = priorityDictionary()   # est.dist. of non-final vert.
		Q[start] = 0
		vOrder = []
		paths = []
		pathlens = []
		
		for v in Q:
			D[v] = Q[v]
			(vx,vy) = setFloat(v)
			
			xIndex = approxWhere(xvalues,vx,10E-3)
			yIndex = approxWhere(yvalues,vy,10E-3)
			x = xIndex[0][0]
			y = yIndex[0][0]
			
			if (endpts.count(v) > 0):
				
				endpts.remove(v)
				Path = shortestPath(D,P,start,v)
				vOrder.append(v)
				paths.append(Path)
				pathlens.append(D[v])
			if  len(endpts) == 0: break
			
			# w is a neighboring vertex
			for w in nWeightPairs[y][x]:
				if str(w[0]) == 'nan':	
					continue
				wStr = setStr((w[0],w[1]))
				vwLength = D[v] + w[2]
				if vwLength < float(edgethreshold):
					if wStr in D:	
						if vwLength < D[wStr]:
							raise ValueError("Dijkstra: found better path to already-final vertex")
				
					elif wStr not in Q or vwLength < Q[wStr]:
						Q[wStr] = vwLength
						P[wStr] = v
			
		return (paths, vOrder,pathlens)		
		
		# End::Dijkstrathreshold
	
# --------------------------------------------------------------------------
	def DijkstraAllPaths(edgethreshold,nWeightPairs,xvalues,\
			yvalues,start,endpts=None):
		"""
		Find shortest paths from the start vertex to all
		vertices nearer than or equal to the end.
		
		Except if path is beyond edgethreshold length.
		Then return string 'Path beyond edge threshold'.

		The input graph G is assumed to have the following
		representation: A vertex can be any object that can
		be used as an index into a dictionary.  G is a
		dictionary, indexed by vertices.  For any vertex v,
		G[v] is itself a dictionary, indexed by the neighbors
		of v.  For any edge v->w, G[v][w] is the length of
		the edge.  This is related to the representation in
		<http://www.python.org/doc/essays/graphs.html>
		where Guido van Rossum suggests representing graphs
		as dictionaries mapping vertices to lists of neighbors,
		however dictionaries of edges have many advantages
		over lists: they can store extra information (here,
		the lengths), they support fast existence tests,
		and they allow easy modification of the graph by edge
		insertion and removal.  Such modifications are not
		needed here but are important in other graph algorithms.
		Since dictionaries obey iterator protocol, a graph
		represented as described here could be handed without
		modification to an algorithm using Guido's representation.

		Of course, G and G[v] need not be Python dict objects;
		they can be any other object that obeys dict protocol,
		for instance a wrapper in which vertices are URLs
		and a call to G[v] loads the web page and finds its links.
		
		The output is a pair (D,P) where D[v] is the distance
		from start to v and P[v] is the predecessor of v along
		the shortest path from s to v.
		
		Dijkstra's algorithm is only guaranteed to work correctly
		when all edge lengths are positive. This code does not
		verify this property for all edges (only the edges seen
		before the end vertex is reached), but will correctly
		compute shortest paths even for some graphs with negative
		edges, and will raise an exception if it discovers that
		a negative edge has caused it to make a mistake.
		"""

		D = {}	# dictionary of final distances
		P = {}	# dictionary of predecessors
		Q = priorityDictionary()   # est.dist. of non-final vert.
		Q[start] = 0
		vOrder = []
		paths = []
		pathlens = []
		xstart = approxWhere(xvalues,float(start.split('_')[0]),10E-3)
		ystart = approxWhere(yvalues,float(start.split('_')[1]),10E-3)
		
		for v in Q:
			D[v] = Q[v]
			#print D[v]		
			
			
			(vx,vy) = setFloat(v)
			
			xIndex = approxWhere(xvalues,vx,10E-3)
			yIndex = approxWhere(yvalues,vy,10E-3)
			x = xIndex[0][0]
			y = yIndex[0][0]
		
			#if x == 51 & y == 17:
			#	print 'distance at endpoint', D[v]
			# w is a neighboring vertex
			for w in nWeightPairs[y][x]:
				if str(w[0]) == 'nan':	
					continue
				wStr = setStr((w[0],w[1]))
				vwLength = D[v] + w[2]
				if vwLength < float(edgethreshold):
					if wStr in D:	
						if vwLength < D[wStr]:
							raise ValueError("Dijkstra: found better path to already-final vertex")
							
					elif wStr not in Q or vwLength < Q[wStr]:
						Q[wStr] = vwLength
						P[wStr] = v
						#if x == 51 & y == 17:
                                                #	print 'distance at endpoint', vwLength
			
		return (D)		
		
		# End::Dijkstrathreshold
	
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
def RasterToNWeights(resgrid):
	'''
	RasterToNWeights()
	Convert raster ascii to Neighbor Weights graph
	'''

	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])

	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)


	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	# Turn no data spot into a really big number - idea is that shortest path would never even look at this value.

	rastervalues[rastervalues == NODATA_value] = np.max(rastervalues) * 1000000
	#rastervalues[rastervalues == NODATA_value] = 1000000

	# Create array
	dt = np.dtype((np.float,(3,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 
	# These common sub-expressions are constants but are 
	# repeated many times in the
	# iteration logic below: isolate these to compute once-only,
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
			
			#Map neighbor weights to each position, setting those positions that don't exist to a null 
			#value, positions are ordered starting from top left corner at zero, and reading from left to 
			#right			
			
			#Position 0: Check on first row, on first column
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
			else:
				zero_xspot = float(xspot-cellsize)
				zero_yspot = float(yspot+cellsize)
				zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
				
			#Position 1: Check on first row
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
			else:
				one_xspot = float(xspot)
				one_yspot = float(yspot+cellsize)
				one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
				
			#Position 2: Check on first row, on last column
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
			else:
				two_xspot = float(xspot+cellsize)
				two_yspot = float(yspot+cellsize)
				two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
			
			#Position 3: Check on first column
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
			else:
				three_xspot = float(xspot-cellsize)
				three_yspot = float(yspot)
				three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
				
			#Position 4: Check on last column
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
			else:
				four_xspot = float(xspot+cellsize)
				four_yspot = float(yspot)
				four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
				
			#Position 5: Check last row, first column 
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None
			
			else:
				five_xspot = float(xspot-cellsize)
				five_yspot = float(yspot-cellsize)
				five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
			
			#Position 6: Check on last row
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
			else:
				six_xspot = float(xspot)
				six_yspot = float(yspot-cellsize)
				six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell	
				
			#Position 7: Check last column, last row
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
			else:
				seven_xspot = float(xspot+cellsize)
				seven_yspot = float(yspot-cellsize)
				seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
			
			# Get neighbors:			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven)	
			
		# end:: irow loop	

	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToNWeights()

# --------------------------------------------------------------------------
def RasterToNWeights2(resgrid):
	'''
	RasterToNWeights()
	Convert raster ascii to Neighbor Weights graph
	'''

	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])

	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)


	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	# Turn no data spot into a really big number - idea is that shortest path would never even look at this value.

	#rastervalues[rastervalues == NODATA_value] = np.max(rastervalues) * 1000000
	#rastervalues[rastervalues == NODATA_value] = 1000000

	# Create array
	dt = np.dtype((np.float,(3,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 
	# These common sub-expressions are constants but are 
	# repeated many times in the
	# iteration logic below: isolate these to compute once-only,
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
			
			#Map neighbor weights to each position, setting those positions that don't exist to a null 
			#value, positions are ordered starting from top left corner at zero, and reading from left to 
			#right			
			
			#Position 0: Check on first row, on first column
			#-----------------------------------------------
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					zero = None
					zero_xspot = None
					zero_yspot = None
				else:
					# Check if zero is nodata, also can't go this way
					if rastervalues[irow-1][icol-1] == -9999:
						zero = None
						zero_xspot = None
						zero_yspot = None
					else:
						zero_xspot = float(xspot-cellsize)
						zero_yspot = float(yspot+cellsize)
						zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
				
			#Position 1: Check on first row
			#------------------------------
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					one = None
					one_xspot = None
					one_yspot = None
				else:
					# Check if one is nodata, also can't go this way
					if rastervalues[irow-1][icol] == -9999:
						one = None
						one_xspot = None
						one_yspot = None
					else:
						one_xspot = float(xspot)
						one_yspot = float(yspot+cellsize)
						one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
				
			#Position 2: Check on first row, on last column
			#----------------------------------------------
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					two = None
					two_xspot = None
					two_yspot = None
				else:
					# Check if two is nodata, also can't go this way
					if rastervalues[irow-1][icol+1] == -9999:
						two_xspot = None
						two_yspot = None
						two = None
					else:
						two_xspot = float(xspot+cellsize)
						two_yspot = float(yspot+cellsize)
						two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
			
			#Position 3: Check on first column
			#---------------------------------
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					three = None
					three_xspot = None
					three_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol-1] == -9999:
						three = None
						three_xspot = None
						three_yspot = None
					else:
						three_xspot = float(xspot-cellsize)
						three_yspot = float(yspot)
						three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
				
			#Position 4: Check on last column
			#--------------------------------
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					four = None
					four_xspot = None
					four_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol+1] == -9999:
						four_xspot = None
						four_yspot = None
						four = None
					else:
						four_xspot = float(xspot+cellsize)
						four_yspot = float(yspot)
						four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
				
			#Position 5: Check last row, first column
			#----------------------------------------
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None
			
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					five = None
					five_xspot = None
					five_yspot = None
				else:
					# Check nodata, also can't go this way
					if rastervalues[irow+1][icol-1] == -9999:
						five_xspot = None
						five_yspot = None
						five = None	
					else:
						five_xspot = float(xspot-cellsize)
						five_yspot = float(yspot-cellsize)
						five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
			
			#Position 6: Check on last row
			#-----------------------------
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					six = None
					six_xspot = None
					six_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol] == -9999:
						six_xspot = None
						six_yspot = None
						six = None
					else:
						six_xspot = float(xspot)
						six_yspot = float(yspot-cellsize)
						six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell	
				
			#Position 7: Check last column, last row
			#---------------------------------------
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					seven = None
					seven_xspot = None
					seven_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol+1] == -9999:
						seven = None
						seven_xspot = None
						seven_yspot = None
					else:
						seven_xspot = float(xspot+cellsize)
						seven_yspot = float(yspot-cellsize)
						seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
			
			# Get neighbors:			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven)	
			
		# end:: irow loop	
	
	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToNWeights()
	

# --------------------------------------------------------------------------
def RasterToBarrierDirection(resgrid,elevgrid,barrgrid,dirtype):
	'''
	RasterToHeightDifference()
	Convert elevation raster ascii to Neighbor height differences
	'''
	
	# Load resistance file
	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])
	
	# Load barrier file
	header_dict_barr, data_list_barr = FileIO.loadFile(barrgrid,header_lines=6)
	
	nrows_barr = int(header_dict_barr['nrows'])
	ncols_barr = int(header_dict_barr['ncols'])
	xllcorner_barr = float(header_dict_barr['xllcorner'])
	yllcorner_barr = float(header_dict_barr['yllcorner'])
	cellsize_barr = float(header_dict_barr['cellsize'])
	NODATA_value_barr = float(header_dict_barr['NODATA_value'])
	
	# Load elevation file
	header_dict_elev, data_list_elev = FileIO.loadFile(elevgrid,header_lines=6)
	
	nrows_elev = int(header_dict_elev['nrows'])
	ncols_elev = int(header_dict_elev['ncols'])
	xllcorner_elev = float(header_dict_elev['xllcorner'])
	yllcorner_elev = float(header_dict_elev['yllcorner'])
	cellsize_elev = float(header_dict_elev['cellsize'])
	NODATA_value_elev = float(header_dict_elev['NODATA_value'])
	
	# Check here to make sure all the same sizes
	if nrows != nrows_barr != nrows_elev and ncols != ncols_barr != ncols_elev:
		print('Elevation and barrier file are not same. (nrows x ncols)')
		sys.exit(-1)
	if xllcorner != xllcorner_barr != xllcorner_elev and yllcorner != yllcorner_barr != yllcorner_elev:
		print('Elevation and barrier file are not same. (xllcorner and yllcorner)')
		sys.exit(-1)
	if cellsize != cellsize_barr != cellsize_elev:
		print('Elevation and barrier file are not same. (cellsize)')
		sys.exit(-1)
		
	# Error check on No data value, can not be -128, this value refers to complete barrier location
	if NODATA_value_barr == -128.0:
		print('No data value can not be -128, this refers to a complete barrier location and cost paths will not traverse this cell.')
	
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)

	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	rastervalues_barr = np.asarray(data_list_barr,dtype = 'float')
	rastervalues_elev = np.asarray(data_list_elev,dtype = 'float')
	
	# Turn no data spot into a really big number - idea is that shortest path would never even look at this value.
	#rastervalues[rastervalues == NODATA_value] = np.max(rastervalues) * 1000000
	#rastervalues[rastervalues == NODATA_value] = 1000000

	
	# Turn no data spot into a really big number - idea is that shortest path would never even look at this value.
	#replaceNA = max(np.max(rastervalues),np.max(rastervalues_barr),np.max(rastervalues_elev)) * 100000
	#replaceNA = None
	#rastervalues[rastervalues == NODATA_value] = replaceNA
	#rastervalues_barr[rastervalues_barr == NODATA_value_barr] = replaceNA
	#rastervalues_elev[rastervalues_elev == NODATA_value_elev] = replaceNA
	
	# Try to leave no data as -9999, then when search around place nan there.
	
	# Note elevation file is now flow accumlation. remove abs. if - then going up if + then going down 
	
	# Note that resistance and barrier file defines extent of streams (overrides elevation file), elvation file can be entire grid not confined to stream extent.
	
	# Create array
	dt = np.dtype((np.float,(3,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 
	# These common sub-expressions are constants but are 
	# repeated many times in the
	# iteration logic below: isolate these to compute once-only,
	# Elevation: Negative is dropping in elevation, positive going up.
	# Flow Accumlation: - is going up, + is dropping down (oposite from elevation)
	# Barrier: Non zero is barrier location.
	# Resistance: Get cost + distance
	# Checks: If elevation < 0 ignore barrier and just calculate cost
	# If elevation >= 0 and barrier > 0, calculate cost + weight from NA replace  
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
			
			#Map neighbor weights to each position, setting those positions that don't exist to a null 
			#value, positions are ordered starting from top left corner at zero, and reading from left to 
			#right
			# If nodata value in neighborhood, place nan			
			
			# Get minimum rastervalues_barr to use when on a barrier replacement
			replacebarrval = np.min(rastervalues_barr.flatten()[np.where(rastervalues_barr.flatten() >= 0)[0]])
			
			# Complete barrier check at this spot
			if rastervalues_barr[irow][icol] == -128.0:
				on_compbarr = True
			else:
				on_compbarr = False
			# Partial barrier check at this spot
			if rastervalues_barr[irow][icol] > 1.0:
				on_partbarr = True
			else:
				on_partbarr = False 
				
			#Position 0: rastervalues_barr[irow-1][icol-1]
			# ---------------------------------------
			# Special case for first row, on first column
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
			else:				
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					zero = None
					zero_xspot = None
					zero_yspot = None
				else:				
					# Check if zero is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol-1] == -9999:
						zero = None
						zero_xspot = None
						zero_yspot = None
					else:					
						# Elevation check
						zero_elev = (rastervalues_elev[irow-1][icol-1]-rastervalues_elev[irow][icol]) 
						# Complete barrier check in zero direction
						if rastervalues_barr[irow-1][icol-1] == -128.0:
							zero_compbarr = True
						else:
							zero_compbarr = False
						# Partial barrier check in zero direction
						if rastervalues_barr[irow-1][icol-1] > 1.0:
							zero_partbarr = True
						else:
							zero_partbarr = False 
											
						#zero_barr = abs(rastervalues_barr[irow-1][icol-1]-rastervalues_barr[irow][icol])
												
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (zero_elev <= 0 and dirtype == 'FlowAcc') or (zero_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									zero = None 
									zero_xspot = None
									zero_yspot = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if zero_compbarr:
										# Can not go this way
										zero = None 
										zero_xspot = None
										zero_yspot = None
									else: 	
										# Can go this way if anything else
										zero = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol-1]) * sqrtOf2xMidCell
										zero_xspot = float(xspot-cellsize)
										zero_yspot = float(yspot+cellsize)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if zero_compbarr or zero_partbarr: 
									zero = (replacebarrval*2) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
								# If not a barrier this way
								else:
									zero = (rastervalues_barr[irow-1][icol-1]+replacebarrval) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
												
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if zero_compbarr or zero_partbarr: 
								# Going up
								if (zero_elev <= 0 and dirtype == 'FlowAcc') or (zero_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if zero_compbarr: # Can't go this way
										zero = None 
										zero_xspot = None
										zero_yspot = None
									elif zero_partbarr: # Can go this way, use barrier costs
										zero = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol-1]) * sqrtOf2xMidCell
										zero_xspot = float(xspot-cellsize)
										zero_yspot = float(yspot+cellsize)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 
									zero = (rastervalues_barr[irow][icol]+replacebarrval) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
							
							else: # No barrier this way
								zero = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol-1]) * sqrtOf2xMidCell
								zero_xspot = float(xspot-cellsize)
								zero_yspot = float(yspot+cellsize)					
					
			#Position 1: rastervalues_barr[irow-1][icol]
			#--------------------------------------
			# Special case for first row
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					one = None
					one_xspot = None
					one_yspot = None
				else:				
					# Check if one is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol] == -9999:
						one = None
						one_xspot = None
						one_yspot = None
					else:						
						# Elevation check
						one_elev = (rastervalues_elev[irow-1][icol]-rastervalues_elev[irow][icol])
						# Complete Barrier check in one direction (up)
						if rastervalues_barr[irow-1][icol] == -128.0:
							one_compbarr = True
						else:
							one_compbarr = False
						# Partial barrier check in zero direction
						if rastervalues_barr[irow-1][icol] > 1.0:
							one_partbarr = True
						else:
							one_partbarr = False
						#one_barr = abs(rastervalues_barr[irow-1][icol]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (one_elev <= 0 and dirtype == 'FlowAcc') or (one_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									one = None
									one_xspot = None
									one_yspot = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if one_compbarr:
										# Can not go this way
										one = None
										one_xspot = None
										one_yspot = None
									else: 	
										# Can go this way if anything else
										one = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol]) * midCell
										one_xspot = float(xspot)
										one_yspot = float(yspot+cellsize)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: # Going down, can go this way
								# If barrier this way too, use replaceval
								if one_compbarr or one_partbarr: # If barrier this way too.
									one = (replacebarrval*2) * midCell 
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
								else:
									one = (replacebarrval+rastervalues_barr[irow-1][icol]) * midCell
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
							
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way						
							if one_compbarr or one_partbarr: # Barrier exists this way
								# Going up
								#---------
								if (one_elev <= 0 and dirtype == 'FlowAcc') or (one_elev >= 0 and dirtype == 'DEM'):  
									# Complete barrier this way
									if one_compbarr:
										one = None
										one_xspot = None
										one_yspot = None
									elif one_partbarr: # Can go this way, use barrier costs
										one = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol]) * midCell
										one_xspot = float(xspot)
										one_yspot = float(yspot+cellsize)	
									else:
										print('Barrier values off.')
										sys.exit(-1)	
								# Going down, can go this way, use replace vals
								else: # Going down, can go this way
									one = (rastervalues_barr[irow][icol]+replacebarrval) * midCell
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
								
							else: # No barrier this way
								one = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol]) * midCell
								one_xspot = float(xspot)
								one_yspot = float(yspot+cellsize)				
				
			#Position 2: rastervalues_elev[irow-1][icol+1]
			#---------------------------------------------
			# Special case for first row, on last column
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					two = None
					two_xspot = None
					two_yspot = None
				else:				
					# Check if two is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol+1] == -9999:
						two_xspot = None
						two_yspot = None
						two = None
					else:
						# Elevation check
						two_elev = (rastervalues_elev[irow-1][icol+1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in two direction
						if rastervalues_barr[irow-1][icol+1] == -128.0:
							two_compbarr = True
						else:
							two_compbarr = False
						# Partial Barrier check in two direction
						if rastervalues_barr[irow-1][icol+1] > 1.0:
							two_partbarr = True
						else:
							two_partbarr = False
						#two_barr = abs(rastervalues_barr[irow-1][icol+1]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							if dirtype == 'FlowAcc':
								# Going up
								#---------
								if (two_elev <= 0 and dirtype == 'FlowAcc') or (two_elev >= 0 and dirtype == 'DEM'): 
									# Comnplete barrier, can't go this way
									if on_compbarr: 
										two_xspot = None
										two_yspot = None
										two = None
									# Can maybe go this way if partial barrier, use partial barrier cost
									elif on_partbarr:
										# If complete barrier this way
										if two_compbarr:
											two_xspot = None
											two_yspot = None
											two = None
										else: 	
											# Can go this way if anything else
											two = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol+1]) * sqrtOf2xMidCell
											two_xspot = float(xspot+cellsize)
											two_yspot = float(yspot+cellsize)
									else:
										print('Barrier raster values not correct.')
										sys.exit(-1)		
								# Going down, can go this way
								#-----------
								else: 
									# If barrier this way too, use replaceval
									if two_compbarr or two_partbarr: # If barrier this way too.
										two = (replacebarrval*2) * sqrtOf2xMidCell
										two_xspot = float(xspot+cellsize)
										two_yspot = float(yspot+cellsize)
									# If not a barrier this way
									else:
										two = (replacebarrval+rastervalues_barr[irow-1][icol+1]) * sqrtOf2xMidCell
										two_xspot = float(xspot+cellsize)
										two_yspot = float(yspot+cellsize)
							
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if two_compbarr or two_partbarr: # Barrier exists this way
								# Going up
								#---------
								if (two_elev <= 0 and dirtype == 'FlowAcc') or (two_elev >= 0 and dirtype == 'DEM'): 
									# Complete barrier this way
									if two_compbarr:
										two_xspot = None
										two_yspot = None
										two = None
									elif two_partbarr: # Can go this way, use barrier costs
										two = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol+1]) * sqrtOf2xMidCell
										two_xspot = float(xspot+cellsize)
										two_yspot = float(yspot+cellsize)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:
									two = (replacebarrval+rastervalues_barr[irow][icol]) * sqrtOf2xMidCell
									two_xspot = float(xspot+cellsize)
									two_yspot = float(yspot+cellsize)
								
							else: # No barrier this way
								two = (rastervalues_barr[irow][icol]+rastervalues_barr[irow-1][icol+1]) * sqrtOf2xMidCell
								two_xspot = float(xspot+cellsize)
								two_yspot = float(yspot+cellsize)
					
			#Position 3: rastervalues_barr[irow][icol-1]
			#-------------------------------------------
			# Special case for first column
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					three = None
					three_xspot = None
					three_yspot = None
				else:
					# Check if three nodata, also can't go this way
					if rastervalues_barr[irow][icol-1] == -9999:
						three = None
						three_xspot = None
						three_yspot = None
					else:
						# Elevation check					
						three_elev = (rastervalues_elev[irow][icol-1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in three direction (left)
						if rastervalues_barr[irow][icol-1] == -128.0:
							three_compbarr = 1
						else:
							three_compbarr = 0
						# Partial Barrier check in three direction (left)
						if rastervalues_barr[irow][icol-1] > 1.0:
							three_partbarr = 1
						else:
							three_partbarr = 0
						#three_barr = abs(rastervalues_barr[irow][icol-1]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (three_elev <= 0 and dirtype == 'FlowAcc') or (three_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									three = None
									three_xspot = None
									three_yspot = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if three_compbarr:
										three = None
										three_xspot = None
										three_yspot = None
									else: 	
										# Can go this way if anything else
										three = (rastervalues_barr[irow][icol]+rastervalues_barr[irow][icol-1]) * midCell
										three_xspot = float(xspot-cellsize)
										three_yspot = float(yspot)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval							
								if three_compbarr or three_partbarr: # If barrier this way too.
									three = (replacebarrval*2) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)

								else:
									three = (rastervalues_barr[irow][icol-1]+replacebarrval) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)			
						
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if three_compbarr or three_partbarr: # Barrier exists this way
								# Going up
								if (three_elev <= 0 and dirtype == 'FlowAcc') or (three_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if three_compbarr:
										three = None
										three_xspot = None
										three_yspot = None
									elif three_partbarr: # Can go this way, use barrier costs
										three = (rastervalues_barr[irow][icol]+rastervalues_barr[irow][icol-1]) * midCell
										three_xspot = float(xspot-cellsize)
										three_yspot = float(yspot)
									else:
										print('Barrier values off.')
										sys.exit(-1)	
								# Going down, can go this way, use replace vals		
								else:		
									three = (rastervalues_barr[irow][icol]+replacebarrval) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)
								
							else: # No barrier this way
								three = (rastervalues_barr[irow][icol]+rastervalues_barr[irow][icol-1]) * midCell
								three_xspot = float(xspot-cellsize)
								three_yspot = float(yspot)						
						
			#Position 4: rastervalues_elev[irow][icol+1]
			#-------------------------------------------
			# Special case for last column
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					four = None
					four_xspot = None
					four_yspot = None
				else:
				
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow][icol+1] == -9999:
						four_xspot = None
						four_yspot = None
						four = None
					else:				
						# Elevation check
						four_elev = (rastervalues_elev[irow][icol+1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in four direction
						if rastervalues_barr[irow][icol+1] == -128.0:
							four_compbarr = True
						else:
							four_compbarr = False
						# Partial Barrier check in four direction
						if rastervalues_barr[irow][icol+1] > 1.0:
							four_partbarr = True
						else:
							four_partbarr = False
						#four_barr = abs(rastervalues_barr[irow][icol+1]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (four_elev <= 0 and dirtype == 'FlowAcc') or (four_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									four_xspot = None
									four_yspot = None
									four = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if four_compbarr:
										four_xspot = None
										four_yspot = None
										four = None
									else: 	
										# Can go this way if anything else
										four = (rastervalues_barr[irow][icol] + rastervalues_barr[irow][icol+1]) * midCell
										four_xspot = float(xspot+cellsize)
										four_yspot = float(yspot)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)		
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if four_compbarr or four_partbarr:
									four = (replacebarrval*2) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)
								# If not a barrier this way
								else:
									four = (replacebarrval + rastervalues_barr[irow][icol+1]) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)			
						
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if four_compbarr or four_partbarr: 
								# Going up
								if (four_elev <= 0 and dirtype == 'FlowAcc') or (four_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if four_compbarr: # Can't go this way
										four_xspot = None
										four_yspot = None
										four = None
									elif four_partbarr: # Can go this way, use barrier costs
										four = (rastervalues_barr[irow][icol] + rastervalues_barr[irow][icol+1]) * midCell
										four_xspot = float(xspot+cellsize)
										four_yspot = float(yspot)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 
									four = (rastervalues_barr[irow][icol] + replacebarrval) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)
								
							else: # No barrier this way
								four = (rastervalues_barr[irow][icol] + rastervalues_barr[irow][icol+1]) * midCell
								four_xspot = float(xspot+cellsize)
								four_yspot = float(yspot)
				
			#Position 5: rastervalues_elev[irow+1][icol-1]
			# --------------------------------------------
			# Special case for last row, first column 
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None			
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					five = None
					five_xspot = None
					five_yspot = None
				else:
				
					# Check nodata, also can't go this way
					if rastervalues_barr[irow+1][icol-1] == -9999:
						five_xspot = None
						five_yspot = None
						five = None	
					else:				
						# Check Elevation
						five_elev = (rastervalues_elev[irow+1][icol-1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in five direction (right)
						if rastervalues_barr[irow+1][icol-1] == -128.0:
							five_compbarr = True
						else:
							five_compbarr = False
						# Partial Barrier check in five direction (right)
						if rastervalues_barr[irow+1][icol-1] > 1.0:
							five_partbarr = True
						else:
							five_partbarr = False
						#five_barr = abs(rastervalues_barr[irow+1][icol-1]-rastervalues_barr[irow][icol])

						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (five_elev <= 0 and dirtype == 'FlowAcc') or (five_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									five_xspot = None
									five_yspot = None
									five = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if five_compbarr:
										five_xspot = None
										five_yspot = None
										five = None
									else: 	
										# Can go this way if anything else
										five = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol-1]) * sqrtOf2xMidCell
										five_xspot = float(xspot-cellsize)
										five_yspot = float(yspot-cellsize)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if five_compbarr or five_partbarr:
									five = (replacebarrval*2) * sqrtOf2xMidCell
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
							
								# If not a barrier this way
								else:
									five = (replacebarrval + rastervalues_barr[irow+1][icol-1]) * sqrtOf2xMidCell	
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if five_compbarr or five_partbarr: 
								# Going up
								if (five_elev <= 0 and dirtype == 'FlowAcc') or (five_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if five_compbarr: # Can't go this way
										five_xspot = None
										five_yspot = None
										five = None	
									elif five_partbarr: # Can go this way, use barrier costs
										five = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol-1]) * sqrtOf2xMidCell
										five_xspot = float(xspot-cellsize)
										five_yspot = float(yspot-cellsize)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:
									five = (rastervalues_barr[irow][icol] + replacebarrval) * sqrtOf2xMidCell	
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
								
							else: # No barrier this way
								five = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol-1]) * sqrtOf2xMidCell
								five_xspot = float(xspot-cellsize)
								five_yspot = float(yspot-cellsize)
						
			#Position 6: rastervalues_elev[irow+1][icol]
			#-------------------------------------------
			# Special case for last row
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					six = None
					six_xspot = None
					six_yspot = None
				else:
				
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow+1][icol] == -9999:
						six_xspot = None
						six_yspot = None
						six = None
					else:
						# Check elevation
						six_elev = (rastervalues_elev[irow+1][icol]-rastervalues_elev[irow][icol])	
						# Complete barrier check in six direction
						if rastervalues_barr[irow+1][icol] == -128.0:
							six_compbarr = True
						else:
							six_compbarr = False
						# Complete barrier check in six direction
						if rastervalues_barr[irow+1][icol] > 1.0:
							six_partbarr = True
						else:
							six_partbarr = False
						#six_barr = (rastervalues_barr[irow+1][icol]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (six_elev <= 0 and dirtype == 'FlowAcc') or (six_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									six_xspot = None
									six_yspot = None
									six = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if six_compbarr:
										six_xspot = None
										six_yspot = None
										six = None
									else: 	
										# Can go this way if anything else
										six = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol]) * midCell
										six_xspot = float(xspot)
										six_yspot = float(yspot-cellsize)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if six_compbarr or six_partbarr:
									six = (repalcebarrva*2) * midCell
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
								# If not a barrier this way
								else:
									six = (replacebarrval + rastervalues_barr[irow+1][icol]) * midCell	
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
											
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if six_compbarr or six_partbarr: 
								# Going up
								if (six_elev <= 0 and dirtype == 'FlowAcc') or (six_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if six_compbarr: # Can't go this way
										six_xspot = None
										six_yspot = None
										six = None
									elif six_partbarr: # Can go this way, use barrier costs
										six = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol]) * midCell
										six_xspot = float(xspot)
										six_yspot = float(yspot-cellsize)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:		
									six = (rastervalues_barr[irow][icol] + replacebarrval) * midCell	
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
								
							else: # No barrier this way
								six = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol]) * midCell
								six_xspot = float(xspot)
								six_yspot = float(yspot-cellsize)
						
			#Position 7: rastervalues_elev[irow+1][icol+1]
			#---------------------------------------------
			# Special case for last column, last row
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					seven = None
					seven_xspot = None
					seven_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow+1][icol+1] == -9999:
						seven = None
						seven_xspot = None
						seven_yspot = None
					else:				
						# Check elevation
						seven_elev = (rastervalues_elev[irow+1][icol+1]-rastervalues_elev[irow][icol])
						# Check complete barrier in seven direction
						if rastervalues_barr[irow+1][icol+1] == -128.0:
							seven_compbarr = True
						else:
							seven_compbarr = False
						# Check complete barrier in seven direction
						if rastervalues_barr[irow+1][icol+1] > 1.0:
							seven_partbarr = True
						else:
							seven_partbarr = False
						#seven_barr = abs(rastervalues_barr[irow+1][icol+1]-rastervalues_barr[irow][icol])
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (seven_elev <= 0 and dirtype == 'FlowAcc') or (seven_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									seven = None
									seven_xspot = None
									seven_yspot = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if seven_compbarr:
										seven = None
										seven_xspot = None
										seven_yspot = None
									else: 	
										# Can go this way if anything else	
										seven = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol+1]) * sqrtOf2xMidCell
										seven_xspot = float(xspot+cellsize)
										seven_yspot = float(yspot-cellsize)
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if seven_compbarr or seven_partbarr: 	
									seven = (replacebarrval*2) * sqrtOf2xMidCell
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)									
								# If not a barrier this way
								else:
									seven = (rastervalues_barr[irow+1][icol+1] + replacebarrval) * sqrtOf2xMidCell	
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)
							
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if seven_compbarr or seven_partbarr: 
								# Going up
								if (seven_elev <= 0 and dirtype == 'FlowAcc') or (seven_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if seven_compbarr: # Can't go this way
										seven = None
										seven_xspot = None
										seven_yspot = None
									elif seven_partbarr: # Can go this way, use barrier costs
										seven = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol+1]) * sqrtOf2xMidCell
										seven_xspot = float(xspot+cellsize)
										seven_yspot = float(yspot-cellsize)
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 	
									seven = (replacebarrval + rastervalues_barr[irow][icol]) * sqrtOf2xMidCell	
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)
								
							else: # No barrier this way
								seven = (rastervalues_barr[irow][icol] + rastervalues_barr[irow+1][icol+1]) * sqrtOf2xMidCell
								seven_xspot = float(xspot+cellsize)
								seven_yspot = float(yspot-cellsize)
				
			# Get neighbors:			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven)	
			
		# end:: irow loop	
	
	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToDirection()
	
# --------------------------------------------------------------------------
def RasterToBarrierDirection_Conductance(resgrid,elevgrid,barrgrid,dirtype):
	'''
	RasterToHeightDifference()
	Convert elevation raster ascii to Neighbor height differences
	'''
	
	# Load resistance file
	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])
	
	# Load barrier file
	header_dict_barr, data_list_barr = FileIO.loadFile(barrgrid,header_lines=6)
	
	nrows_barr = int(header_dict_barr['nrows'])
	ncols_barr = int(header_dict_barr['ncols'])
	xllcorner_barr = float(header_dict_barr['xllcorner'])
	yllcorner_barr = float(header_dict_barr['yllcorner'])
	cellsize_barr = float(header_dict_barr['cellsize'])
	NODATA_value_barr = float(header_dict_barr['NODATA_value'])
	
	# Load elevation file
	header_dict_elev, data_list_elev = FileIO.loadFile(elevgrid,header_lines=6)
	
	nrows_elev = int(header_dict_elev['nrows'])
	ncols_elev = int(header_dict_elev['ncols'])
	xllcorner_elev = float(header_dict_elev['xllcorner'])
	yllcorner_elev = float(header_dict_elev['yllcorner'])
	cellsize_elev = float(header_dict_elev['cellsize'])
	NODATA_value_elev = float(header_dict_elev['NODATA_value'])
	
	# Check here to make sure all the same sizes
	if nrows != nrows_barr != nrows_elev and ncols != ncols_barr != ncols_elev:
		print('Elevation and barrier file are not same. (nrows x ncols)')
		sys.exit(-1)
	if xllcorner != xllcorner_barr != xllcorner_elev and yllcorner != yllcorner_barr != yllcorner_elev:
		print('Elevation and barrier file are not same. (xllcorner and yllcorner)')
		sys.exit(-1)
	if cellsize != cellsize_barr != cellsize_elev:
		print('Elevation and barrier file are not same. (cellsize)')
		sys.exit(-1)
		
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)

	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	rastervalues_barr = np.asarray(data_list_barr,dtype = 'float')
	rastervalues_elev = np.asarray(data_list_elev,dtype = 'float')
			
	# Try to leave no data as -9999, then when search around place nan there.
	
	# Note elevation file is now flow accumlation. remove abs. if - then going up if + then going down 
	
	# Note that resistance and barrier file defines extent of streams (overrides elevation file), elvation file can be entire grid not confined to stream extent.
	
	# Note no longer using resgrid:
	
	# Still need resistance for direction, if using conductance then 1/conductance to get resistance. if 0 conducation assume complete barrier
	
	# Create array: x,y,resistance,conductance
	dt = np.dtype((np.float,(4,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 
	# These common sub-expressions are constants but are 
	# repeated many times in the
	# iteration logic below: isolate these to compute once-only,
	# Elevation: Negative is dropping in elevation, positive going up.
	# Flow Accumlation: - is going up, + is dropping down (oposite from elevation)
	# Barrier: Non zero is barrier location.
	# Resistance: Get cost + distance
	# Checks: If elevation < 0 ignore barrier and just calculate cost
	# If elevation >= 0 and barrier > 0, calculate cost + weight from NA replace  
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
			
			#Map neighbor weights to each position, setting those positions that don't exist to a null 
			#value, positions are ordered starting from top left corner at zero, and reading from left to 
			#right
			# If nodata value in neighborhood, place nan			
			
			# Get minimum rastervalues_barr to use when on a barrier replacement
			#replacebarrval = np.min(rastervalues_barr.flatten()[np.where(rastervalues_barr.flatten() >= 0)[0]])
			replacebarrval = 1.0 # Assume that it is 1!
			
			# Complete barrier check at this spot
			if rastervalues_barr[irow][icol] == 0.0:
				on_compbarr = True
			else:
				on_compbarr = False
			# Partial barrier check at this spot
			if 0.0 < rastervalues_barr[irow][icol] < 1.0:
				on_partbarr = True
			else:
				on_partbarr = False 
			
			#Position 0: rastervalues_barr[irow-1][icol-1]
			# ---------------------------------------
			# Special case for first row, on first column
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
				zero_cond = None
			else:				
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					zero = None
					zero_xspot = None
					zero_yspot = None
					zero_cond = None
				else:				
					# Check if zero is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol-1] == -9999:
						zero = None
						zero_xspot = None
						zero_yspot = None
						zero_cond = None
					else:					
						# Elevation check
						zero_elev = (rastervalues_elev[irow-1][icol-1]-rastervalues_elev[irow][icol]) 
						# Complete barrier check in zero direction
						if rastervalues_barr[irow-1][icol-1] == 0.0:
							zero_compbarr = True
						else:
							zero_compbarr = False
						# Partial barrier check in zero direction
						if  0.0 < rastervalues_barr[irow-1][icol-1] < 1.0:
							zero_partbarr = True
						else:
							zero_partbarr = False 
											
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (zero_elev <= 0 and dirtype == 'FlowAcc') or (zero_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									zero = None 
									zero_xspot = None
									zero_yspot = None
									zero_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if zero_compbarr:
										# Can not go this way
										zero = None 
										zero_xspot = None
										zero_yspot = None
										zero_cond = None
									else: 	
										# Can go this way if anything else
										zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
										zero_xspot = float(xspot-cellsize)
										zero_yspot = float(yspot+cellsize)
										zero_cond = rastervalues_barr[irow-1][icol-1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if zero_compbarr or zero_partbarr: 
									zero = (replacebarrval*2) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
									zero_cond = replacebarrval
								# If not a barrier this way
								else:
									zero = (rastervalues[irow-1][icol-1]+replacebarrval) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
									zero_cond = replacebarrval
												
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if zero_compbarr or zero_partbarr: 
								# Going up
								if (zero_elev <= 0 and dirtype == 'FlowAcc') or (zero_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if zero_compbarr: # Can't go this way
										zero = None 
										zero_xspot = None
										zero_yspot = None
										zero_cond = None
									elif zero_partbarr: # Can go this way, use barrier costs
										zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
										zero_xspot = float(xspot-cellsize)
										zero_yspot = float(yspot+cellsize)
										zero_cond = rastervalues_barr[irow-1][icol-1]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 
									zero = (rastervalues[irow][icol]+replacebarrval) * sqrtOf2xMidCell
									zero_xspot = float(xspot-cellsize)
									zero_yspot = float(yspot+cellsize)
									zero_cond = replacebarrval
							
							else: # No barrier this way
								zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
								zero_xspot = float(xspot-cellsize)
								zero_yspot = float(yspot+cellsize)
								zero_cond = rastervalues_barr[irow-1][icol-1]								
					
			#Position 1: rastervalues_barr[irow-1][icol]
			#--------------------------------------
			# Special case for first row
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
				one_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					one = None
					one_xspot = None
					one_yspot = None
					one_cond = None
				else:				
					# Check if one is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol] == -9999:
						one = None
						one_xspot = None
						one_yspot = None
						one_cond = None
					else:						
						# Elevation check
						one_elev = (rastervalues_elev[irow-1][icol]-rastervalues_elev[irow][icol])
						# Complete Barrier check in one direction (up)
						if rastervalues_barr[irow-1][icol] == 0.0:
							one_compbarr = True
						else:
							one_compbarr = False
						# Partial barrier check in zero direction
						if 0.0 < rastervalues_barr[irow-1][icol] < 1.0:
							one_partbarr = True
						else:
							one_partbarr = False
												
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (one_elev <= 0 and dirtype == 'FlowAcc') or (one_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									one = None
									one_xspot = None
									one_yspot = None
									one_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if one_compbarr:
										# Can not go this way
										one = None
										one_xspot = None
										one_yspot = None
										one_cond = None
									else: 	
										# Can go this way if anything else
										one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
										one_xspot = float(xspot)
										one_yspot = float(yspot+cellsize)
										one_cond = rastervalues_barr[irow-1][icol]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: # Going down, can go this way
								# If barrier this way too, use replaceval
								if one_compbarr or one_partbarr: # If barrier this way too.
									one = (replacebarrval*2) * midCell 
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
									one_cond = replacebarrval
								else:
									one = (replacebarrval+rastervalues[irow-1][icol]) * midCell
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
									one_cond = replacebarrval
							
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way						
							if one_compbarr or one_partbarr: # Barrier exists this way
								# Going up
								#---------
								if (one_elev <= 0 and dirtype == 'FlowAcc') or (one_elev >= 0 and dirtype == 'DEM'):  
									# Complete barrier this way
									if one_compbarr:
										one = None
										one_xspot = None
										one_yspot = None
									elif one_partbarr: # Can go this way, use barrier costs
										one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
										one_xspot = float(xspot)
										one_yspot = float(yspot+cellsize)
										one_cond = rastervalues_barr[irow-1][icol]										
									else:
										print('Barrier values off.')
										sys.exit(-1)	
								# Going down, can go this way, use replace vals
								else: # Going down, can go this way
									one = (rastervalues[irow][icol]+replacebarrval) * midCell
									one_xspot = float(xspot)
									one_yspot = float(yspot+cellsize)
									one_cond = replacebarrval
								
							else: # No barrier this way
								one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
								one_xspot = float(xspot)
								one_yspot = float(yspot+cellsize)
								one_cond = rastervalues_barr[irow-1][icol]								
				
			#Position 2: rastervalues_elev[irow-1][icol+1]
			#---------------------------------------------
			# Special case for first row, on last column
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
				two_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					two = None
					two_xspot = None
					two_yspot = None
					two_cond = None
				else:				
					# Check if two is nodata, also can't go this way
					if rastervalues_barr[irow-1][icol+1] == -9999:
						two_xspot = None
						two_yspot = None
						two = None
						two_cond = None
					else:
						# Elevation check
						two_elev = (rastervalues_elev[irow-1][icol+1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in two direction
						if rastervalues_barr[irow-1][icol+1] == 0.0:
							two_compbarr = True
							
						else:
							two_compbarr = False
						# Partial Barrier check in two direction
						if 0.0 < rastervalues_barr[irow-1][icol+1] < 1.0:
							two_partbarr = True
						else:
							two_partbarr = False
						#two_barr = abs(rastervalues_barr[irow-1][icol+1]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (two_elev <= 0 and dirtype == 'FlowAcc') or (two_elev >= 0 and dirtype == 'DEM'): 
								# Comnplete barrier, can't go this way
								if on_compbarr: 
									two_xspot = None
									two_yspot = None
									two = None
									two_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if two_compbarr:
										two_xspot = None
										two_yspot = None
										two = None
										two_cond = None
									else: 	
										# Can go this way if anything else
										two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
										two_xspot = float(xspot+cellsize)
										two_yspot = float(yspot+cellsize)
										two_cond = rastervalues_barr[irow-1][icol+1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)		
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if two_compbarr or two_partbarr: # If barrier this way too.
									two = (replacebarrval*2) * sqrtOf2xMidCell
									two_xspot = float(xspot+cellsize)
									two_yspot = float(yspot+cellsize)
									two_cond = replacebarrval
								# If not a barrier this way
								else:
									two = (replacebarrval+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
									two_xspot = float(xspot+cellsize)
									two_yspot = float(yspot+cellsize)
									two_cond = replacebarrval
						
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if two_compbarr or two_partbarr: # Barrier exists this way
								# Going up
								#---------
								if (two_elev <= 0 and dirtype == 'FlowAcc') or (two_elev >= 0 and dirtype == 'DEM'): 
									# Complete barrier this way
									if two_compbarr:
										two_xspot = None
										two_yspot = None
										two = None
										two_cond = None
									elif two_partbarr: # Can go this way, use barrier costs
										two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
										two_xspot = float(xspot+cellsize)
										two_yspot = float(yspot+cellsize)
										two_cond = rastervalues_barr[irow-1][icol+1]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:
									two = (replacebarrval+rastervalues[irow][icol]) * sqrtOf2xMidCell
									two_xspot = float(xspot+cellsize)
									two_yspot = float(yspot+cellsize)
									two_cond = replacebarrval
								
							else: # No barrier this way
								two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
								two_xspot = float(xspot+cellsize)
								two_yspot = float(yspot+cellsize)
								two_cond = rastervalues_barr[irow-1][icol+1]
					
			#Position 3: rastervalues_barr[irow][icol-1]
			#-------------------------------------------
			# Special case for first column
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
				three_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					three = None
					three_xspot = None
					three_yspot = None
					three_cond = None
				else:
					# Check if three nodata, also can't go this way
					if rastervalues_barr[irow][icol-1] == -9999:
						three = None
						three_xspot = None
						three_yspot = None
						three_cond = None
					else:
						# Elevation check					
						three_elev = (rastervalues_elev[irow][icol-1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in three direction (left)
						if rastervalues_barr[irow][icol-1] == 0.0:
							three_compbarr = 1
							
						else:
							three_compbarr = 0
						# Partial Barrier check in three direction (left)
						if 0.0 < rastervalues_barr[irow][icol-1] < 1.0:
							three_partbarr = 1
						else:
							three_partbarr = 0
												
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (three_elev <= 0 and dirtype == 'FlowAcc') or (three_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									three = None
									three_xspot = None
									three_yspot = None
									three_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if three_compbarr:
										three = None
										three_xspot = None
										three_yspot = None
										three_cond = None
									else: 	
										# Can go this way if anything else
										three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
										three_xspot = float(xspot-cellsize)
										three_yspot = float(yspot)
										three_cond = rastervalues_barr[irow][icol-1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval							
								if three_compbarr or three_partbarr: # If barrier this way too.
									three = (replacebarrval*2) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)
									three_cond = replacebarrval

								else:
									three = (rastervalues[irow][icol-1]+replacebarrval) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)	
									three_cond = replacebarrval
						
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if three_compbarr or three_partbarr: # Barrier exists this way
								# Going up
								if (three_elev <= 0 and dirtype == 'FlowAcc') or (three_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if three_compbarr:
										three = None
										three_xspot = None
										three_yspot = None
										three_cond = None
									elif three_partbarr: # Can go this way, use barrier costs
										three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
										three_xspot = float(xspot-cellsize)
										three_yspot = float(yspot)
										three_cond = rastervalues_barr[irow][icol-1]
									else:
										print('Barrier values off.')
										sys.exit(-1)	
								# Going down, can go this way, use replace vals		
								else:		
									three = (rastervalues[irow][icol]+replacebarrval) * midCell
									three_xspot = float(xspot-cellsize)
									three_yspot = float(yspot)
									three_cond = replacebarrval
								
							else: # No barrier this way
								three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
								three_xspot = float(xspot-cellsize)
								three_yspot = float(yspot)
								three_cond = rastervalues_barr[irow][icol-1]								
						
			#Position 4: rastervalues_elev[irow][icol+1]
			#-------------------------------------------
			# Special case for last column
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
				four_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					four = None
					four_xspot = None
					four_yspot = None
					four_cond = None
				else:
				
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow][icol+1] == -9999:
						four_xspot = None
						four_yspot = None
						four = None
						four_cond = None
					else:				
						# Elevation check
						four_elev = (rastervalues_elev[irow][icol+1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in four direction
						if rastervalues_barr[irow][icol+1] == 0.0:
							four_compbarr = True
							
						else:
							four_compbarr = False
						# Partial Barrier check in four direction
						if 0.0 < rastervalues_barr[irow][icol+1] < 1.0:
							four_partbarr = True
						else:
							four_partbarr = False
												
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (four_elev <= 0 and dirtype == 'FlowAcc') or (four_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									four_xspot = None
									four_yspot = None
									four = None
									four_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if four_compbarr:
										four_xspot = None
										four_yspot = None
										four = None
										four_cond = None
									else: 	
										# Can go this way if anything else
										four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
										four_xspot = float(xspot+cellsize)
										four_yspot = float(yspot)
										four_cond = rastervalues_barr[irow][icol+1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)		
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if four_compbarr or four_partbarr:
									four = (replacebarrval*2) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)
									four_cond = replacebarrval
								# If not a barrier this way
								else:
									four = (replacebarrval + rastervalues[irow][icol+1]) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)	
									four_cond = replacebarrval
						
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if four_compbarr or four_partbarr: 
								# Going up
								if (four_elev <= 0 and dirtype == 'FlowAcc') or (four_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if four_compbarr: # Can't go this way
										four_xspot = None
										four_yspot = None
										four = None
										four_cond = None
									elif four_partbarr: # Can go this way, use barrier costs
										four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
										four_xspot = float(xspot+cellsize)
										four_yspot = float(yspot)
										four_cond = rastervalues_barr[irow][icol+1]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 
									four = (rastervalues[irow][icol] + replacebarrval) * midCell
									four_xspot = float(xspot+cellsize)
									four_yspot = float(yspot)
									four_cond = replacebarrval
								
							else: # No barrier this way
								four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
								four_xspot = float(xspot+cellsize)
								four_yspot = float(yspot)
								four_cond = rastervalues_barr[irow][icol+1]
				
			#Position 5: rastervalues_elev[irow+1][icol-1]
			# --------------------------------------------
			# Special case for last row, first column 
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None	
				five_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					five = None
					five_xspot = None
					five_yspot = None
					five_cond = None
				else:
				
					# Check nodata, also can't go this way
					if rastervalues_barr[irow+1][icol-1] == -9999:
						five_xspot = None
						five_yspot = None
						five = None	
						five_cond = None
					else:				
						# Check Elevation
						five_elev = (rastervalues_elev[irow+1][icol-1]-rastervalues_elev[irow][icol])
						# Complete Barrier check in five direction 
						if rastervalues_barr[irow+1][icol-1] == 0.0:
							five_compbarr = True
							
						else:
							five_compbarr = False
						# Partial Barrier check in five direction (right)
						if 0.0 < rastervalues_barr[irow+1][icol-1] < 1.0:
							five_partbarr = True
						else:
							five_partbarr = False
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (five_elev <= 0 and dirtype == 'FlowAcc') or (five_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr:
									five_xspot = None
									five_yspot = None
									five = None
									five_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if five_compbarr:
										five_xspot = None
										five_yspot = None
										five = None
										five_cond = None
									else: 	
										# Can go this way if anything else
										five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
										five_xspot = float(xspot-cellsize)
										five_yspot = float(yspot-cellsize)
										five_cond = rastervalues_barr[irow+1][icol-1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if five_compbarr or five_partbarr:
									five = (replacebarrval*2) * sqrtOf2xMidCell
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
									five_cond = replacebarrval
							
								# If not a barrier this way
								else:
									five = (replacebarrval + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell	
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
									five_cond = replacebarrval
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if five_compbarr or five_partbarr: 
								# Going up
								if (five_elev <= 0 and dirtype == 'FlowAcc') or (five_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if five_compbarr: # Can't go this way
										five_xspot = None
										five_yspot = None
										five = None	
										five_cond = None
									elif five_partbarr: # Can go this way, use barrier costs
										five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
										five_xspot = float(xspot-cellsize)
										five_yspot = float(yspot-cellsize)
										five_cond = rastervalues_barr[irow+1][icol-1]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:
									five = (rastervalues[irow][icol] + replacebarrval) * sqrtOf2xMidCell	
									five_xspot = float(xspot-cellsize)
									five_yspot = float(yspot-cellsize)
									five_cond = replacebarrval
								
							else: # No barrier this way
								five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
								five_xspot = float(xspot-cellsize)
								five_yspot = float(yspot-cellsize)
								five_cond = rastervalues_barr[irow+1][icol-1]
						
			#Position 6: rastervalues_elev[irow+1][icol]
			#-------------------------------------------
			# Special case for last row
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
				six_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					six = None
					six_xspot = None
					six_yspot = None
					six_cond = None
				else:
				
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow+1][icol] == -9999:
						six_xspot = None
						six_yspot = None
						six = None
						six_cond = None
					else:
						# Check elevation
						six_elev = (rastervalues_elev[irow+1][icol]-rastervalues_elev[irow][icol])	
						# Complete barrier check in six direction
						if rastervalues_barr[irow+1][icol] == 0.0:
							six_compbarr = True
							
						else:
							six_compbarr = False
						# Complete barrier check in six direction
						if 0.0 < rastervalues_barr[irow+1][icol] < 1.0:
							six_partbarr = True
						else:
							six_partbarr = False
						#six_barr = (rastervalues_barr[irow+1][icol]-rastervalues_barr[irow][icol])
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (six_elev <= 0 and dirtype == 'FlowAcc') or (six_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									six_xspot = None
									six_yspot = None
									six = None
									six_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if six_compbarr:
										six_xspot = None
										six_yspot = None
										six = None
										six_cond = None
									else: 	
										# Can go this way if anything else
										six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell
										six_xspot = float(xspot)
										six_yspot = float(yspot-cellsize)
										six_cond = rastervalues_barr[irow+1][icol]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if six_compbarr or six_partbarr:
									six = (repalcebarrva*2) * midCell
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
									six_cond = replacebarrval
								# If not a barrier this way
								else:
									six = (replacebarrval + rastervalues[irow+1][icol]) * midCell	
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
									six_cond = replacebarrval
											
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if six_compbarr or six_partbarr: 
								# Going up
								if (six_elev <= 0 and dirtype == 'FlowAcc') or (six_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if six_compbarr: # Can't go this way
										six_xspot = None
										six_yspot = None
										six = None
										six_cond = None
									elif six_partbarr: # Can go this way, use barrier costs
										six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell
										six_xspot = float(xspot)
										six_yspot = float(yspot-cellsize)
										six_cond = rastervalues_barr[irow+1][icol]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else:		
									six = (rastervalues[irow][icol] + replacebarrval) * midCell	
									six_xspot = float(xspot)
									six_yspot = float(yspot-cellsize)
									six_cond = replacebarrval
								
							else: # No barrier this way
								six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell
								six_xspot = float(xspot)
								six_yspot = float(yspot-cellsize)
								six_cond = rastervalues_barr[irow+1][icol]
						
			#Position 7: rastervalues_elev[irow+1][icol+1]
			#---------------------------------------------
			# Special case for last column, last row
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
				seven_cond = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues_barr[irow][icol] == -9999:
					seven = None
					seven_xspot = None
					seven_yspot = None
					seven_cond = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues_barr[irow+1][icol+1] == -9999:
						seven = None
						seven_xspot = None
						seven_yspot = None
						seven_cond = None
					else:				
						# Check elevation
						seven_elev = (rastervalues_elev[irow+1][icol+1]-rastervalues_elev[irow][icol])
						# Check complete barrier in seven direction
						if rastervalues_barr[irow+1][icol+1] == 0.0:
							seven_compbarr = True
							
						else:
							seven_compbarr = False
						# Check complete barrier in seven direction
						if 0.0 < rastervalues_barr[irow+1][icol+1] < 1.0:
							seven_partbarr = True
						else:
							seven_partbarr = False
						
						# Check if on barrier
						#--------------------
						if on_compbarr or on_partbarr:
							# Going up
							#---------
							if (seven_elev <= 0 and dirtype == 'FlowAcc') or (seven_elev >= 0 and dirtype == 'DEM'): 
								# Can't go this way if on complete barrier
								if on_compbarr: 
									seven = None
									seven_xspot = None
									seven_yspot = None
									seven_cond = None
								# Can maybe go this way if partial barrier, use partial barrier cost
								elif on_partbarr:
									# If complete barrier this way
									if seven_compbarr:
										seven = None
										seven_xspot = None
										seven_yspot = None
										seven_cond = None
									else: 	
										# Can go this way if anything else	
										seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
										seven_xspot = float(xspot+cellsize)
										seven_yspot = float(yspot-cellsize)
										seven_cond = rastervalues_barr[irow+1][icol+1]
								else:
									print('Barrier raster values not correct.')
									sys.exit(-1)
							# Going down, can go this way
							#-----------
							else: 
								# If barrier this way too, use replaceval
								if seven_compbarr or seven_partbarr: 	
									seven = (replacebarrval*2) * sqrtOf2xMidCell
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)
									seven_cond = replacebarrval									
								# If not a barrier this way
								else:
									seven = (rastervalues[irow+1][icol+1] + replacebarrval) * sqrtOf2xMidCell	
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)
									seven_cond = replacebarrval
							
						# Not on a barrier
						#-----------------
						else: 							
							# Barrier exists this way
							if seven_compbarr or seven_partbarr: 
								# Going up
								if (seven_elev <= 0 and dirtype == 'FlowAcc') or (seven_elev >= 0 and dirtype == 'DEM'):
									# Complete barrier this way
									if seven_compbarr: # Can't go this way
										seven = None
										seven_xspot = None
										seven_yspot = None
										seven_cond = None
									elif seven_partbarr: # Can go this way, use barrier costs
										seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
										seven_xspot = float(xspot+cellsize)
										seven_yspot = float(yspot-cellsize)
										seven_cond = rastervalues_barr[irow+1][icol+1]
									else:
										print('Barrier values off.')
										sys.exit(-1)
								# Going down, can go this way, use replace vals		
								else: 	
									seven = (replacebarrval + rastervalues[irow][icol]) * sqrtOf2xMidCell	
									seven_xspot = float(xspot+cellsize)
									seven_yspot = float(yspot-cellsize)
									seven_cond = replacebarrval
								
							else: # No barrier this way
								seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
								seven_xspot = float(xspot+cellsize)
								seven_yspot = float(yspot-cellsize)
								seven_cond = rastervalues_barr[irow+1][icol+1]
			if (zero != None and zero < 0) or (one != None and one < 0) or (two != None and two <0) or (three !=None and three <0) or (four != None and four <0) or (five !=None and five <0) or (six != None and six <0) or (seven != None and seven <0):
				print(('Resgrid or barrgrid or facgrid not matching and negative value found around, X='+str(xpot)+' Y='+str(yspot)))
				sys.exit(-1)
			# Get neighbors:			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero,zero_cond)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one,one_cond)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two,two_cond)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three,three_cond)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four,four_cond)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five,five_cond)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six,six_cond)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven,seven_cond)	
			
		# end:: irow loop	
	
	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToDirection()

# --------------------------------------------------------------------------
def RasterToWindSpeedDirection(resgrid,elevgrid,barrgrid,maxres,minres):
	'''
	RasterToWindSpeedDifference()
	Convert U and V grids to resultant wind speed in Neighbor. NWeights  
	then 
	'''
	
	# Load resistance file
	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])
	
	# Load barrier file - this is u direction, rename so not confused.
	header_dict_u, data_list_u = FileIO.loadFile(barrgrid,header_lines=6)
	
	nrows_u = int(header_dict_u['nrows'])
	ncols_u = int(header_dict_u['ncols'])
	xllcorner_u = float(header_dict_u['xllcorner'])
	yllcorner_u = float(header_dict_u['yllcorner'])
	cellsize_u = float(header_dict_u['cellsize'])
	NODATA_value_u = float(header_dict_u['NODATA_value'])
	
	# Load elevation file - this is v direction, rename so not confused.
	header_dict_v, data_list_v = FileIO.loadFile(elevgrid,header_lines=6)
	
	nrows_v = int(header_dict_v['nrows'])
	ncols_v = int(header_dict_v['ncols'])
	xllcorner_v = float(header_dict_v['xllcorner'])
	yllcorner_v = float(header_dict_v['yllcorner'])
	cellsize_v = float(header_dict_v['cellsize'])
	NODATA_value_v = float(header_dict_v['NODATA_value'])
	
	# Check here to make sure all the same sizes
	if nrows != nrows_u != nrows_v and ncols != ncols_u != ncols_v:
		print('U and V file are not same. (nrows x ncols)')
		sys.exit(-1)
	if xllcorner != xllcorner_u != xllcorner_v and yllcorner != yllcorner_u != yllcorner_v:
		print('U and V file are not same. (xllcorner and yllcorner)')
		sys.exit(-1)
	if cellsize != cellsize_u != cellsize_v:
		print('U and V file are not same. (cellsize)')
		sys.exit(-1)
		
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)

	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	rastervalues_u = np.asarray(data_list_u,dtype = 'float')
	rastervalues_v = np.asarray(data_list_v,dtype = 'float')
		
	# Leave no data as -9999, then when search around place nan there.
		
	# Note that resistance (resgrid) defines extent of streams (overrides elevation file), elvation file can be entire grid not confined to stream extent. Leave this for now, may not use resgrid. Maybe include resgrid 
	
	# Sum up the resgrid cost and wind speed. Wind speed direction is largest value. 
	
	# Create array: x,y,resgrid cost, wind speed, 1 - wind direction probability,
	# x,y,Euclidean cost, wind speed - at end will calculate wind speed to resistance
	dt = np.dtype((np.float,(4,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 
	# These common sub-expressions are constants but are 
	# repeated many times in the
	# iteration logic below: isolate these to compute once-only,
	# Resistance: Get cost + distance
	# Wind: Get resultant speed  
	# phi: 0 - 315, 1 - 0, 2 - 45, 3 - 270, 4 - 90, 5 - 225, 6 - 180, 7 - 135
	# Resultant_x = U + ud, where ud = cos(phi)
	# Resultant_y = V + vd, where vd = sin(phi)
	# Speed = sqrt(resultant_x**2 + resultant_y**2)
	#Map neighbor weights to each position, setting those positions that don't exist to a null 
	#value, positions are ordered starting from top left corner at zero, and reading from left to 
	#right
	# If nodata value in neighborhood, place nan	
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	
	# The index and order of the wind direction angles
	phi_ind = [0,1,2,3,4,5,6,7]
	#phi_arr = [315,0,45,270,90,225,180,135] # meterological wind direction
	phi_arr = [135,90,45,180,0,225,270,315] # math direction
	phi_arr = np.asarray(phi_arr) * math.pi / 180.
	r2d = 45./math.atan(1.0)
	#order_cost = np.asarray([1,2,3,4,5,6,7,np.nan])
	order_cost = np.asarray([2,3,4,5,6,7,8,9])
	
	# Begin irow loop
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		# Begin icol loop
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
			
			# Get nearest direction neighbor
			if rastervalues_v[irow][icol] != -9999:
				wind_dir = np.arctan2(rastervalues_v[irow][icol],rastervalues_u[irow][icol]) * 180 / np.pi
				if wind_dir < 0: # correct for atan2 values that are negative
					wind_dir = wind_dir + 360
				takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
				closest_dir = takeClosest(wind_dir,phi_arr)
				closest_ind = np.where(np.asarray(phi_arr) == closest_dir)[0][0]
			else:
				closest_ind = -9999
			
			# Convert neighbors to probabilities
			ccw = wind_dir + 180 # 
			cc = wind_dir - 180
			#ccw_dir = (-1./180) * (phi_arr - wind_dir) + 1
			#cc_dir = (1./180) * (phi_arr - wind_dir) + 1
			wind_prob = []
			for iwind in range(8):
				phi = phi_arr[iwind] # This neighbor
				# For wind_dir 0 to 180, then ccw is less than 360
				if wind_dir <= 180 and wind_dir >= 0:
					# Then if phi falls in the ccw angles
					if phi >= wind_dir and phi <= ccw:
						# Get probability
						wind_prob.append((-1./180) * (phi - wind_dir) + 1)
					else: # Phi not in ccw direction, get corrsponding angle
						alpha = phi - wind_dir
						phi_adj = wind_dir - alpha
						if phi_adj < 0:
							phi_adj = phi_adj +360
						wind_prob.append((-1./180) * (phi_adj - wind_dir) + 1)
						
				# wind_dir is 180 to 360, use cc direction	
				else: 
					# Then if phi falls in the cc angles
					if phi < wind_dir and phi > cc:
						# Get probability
						wind_prob.append((1./180) * (phi - wind_dir) + 1)
					else: # phi not in cc directin, get corresponding angle
						alpha = phi - wind_dir
						phi_adj = wind_dir - alpha
						if phi_adj > 360:
							phi_adj = phi_adj -360
						wind_prob.append((1./180) * (phi_adj - wind_dir) + 1)
			
			# Take 1 - wind prob to use for shortest path
			wind_prob = 1. - np.asarray(wind_prob)
			
			# Get nearest neighbors and order linear with cost
			diff = np.abs(wind_dir - np.asarray(phi_arr))
			sort_index = np.argsort(diff)
			
			# Direction and speed of this wind vector
			md = math.atan2(rastervalues_v[irow][icol],rastervalues_u[irow][icol]) * r2d			
			#wwd = math.atan2(rastervalues_u[irow][icol],rastervalues_v[irow][icol]) * r2d +180
			speed = np.sqrt(rastervalues_u[irow][icol]**2+rastervalues_v[irow][icol]**2)
						
			#Position 0: rastervalues_u[irow-1][icol-1]
			# ---------------------------------------
			phi = phi_arr[0]
			phi_ind = 0
			# Special case for first row, on first column
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
				zero_wind = None
				zero_dir = None
			else:				
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					zero = None
					zero_xspot = None
					zero_yspot = None
					zero_wind = None
					zero_dir = None
				else:				
					# Check if zero is nodata, also can't go this way
					if rastervalues[irow-1][icol-1] == -9999:
						zero = None
						zero_xspot = None
						zero_yspot = None
						zero_wind = None
						zero_dir = None
					else:
						
						# Get wind speed resultant in zero direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						#math.atan2(R_x,R_y) * r2d +180 from direction
						zero_xspot = float(xspot-cellsize)
						zero_yspot = float(yspot+cellsize)
						zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]) * sqrtOf2xMidCell
						zero_wind = Speed
						#zero_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							#zero_dir = cost* sqrtOf2xMidCell
							zero_dir = cost
						else:
							zero = None
							zero_xspot = None
							zero_yspot = None
							zero_wind = None
							zero_dir = None
						
			#------------------------------
			#Position 1: rastervalues_u[irow-1][icol]
			#------------------------------
			phi = 90.
			phi = phi_arr[1]
			phi_ind = 1
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
				one_wind = None
				one_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					one = None
					one_xspot = None
					one_yspot = None
					one_wind = None
					one_dir = None
				else:
					# Check if one is nodata, also can't go this way
					if rastervalues[irow-1][icol] == -9999:
						one = None
						one_xspot = None
						one_yspot = None
						one_wind = None
						one_dir = None
					else:
						# Get wind speed resultant in one direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						one_xspot = float(xspot)
						one_yspot = float(yspot+cellsize)
						one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]) * midCell
						one_wind = Speed
						one_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							one_dir = cost
						else:
							one = None
							one_xspot = None
							one_yspot = None
							one_wind = None
							one_dir = None
						
			# ---------------------------------------------
			#Position 2: Check on first row, on last column [irow-1][icol+1]
			#----------------------------------------------
			phi = 45.
			phi = phi_arr[2]
			phi_ind = 2
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
				two_wind = None
				two_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					two = None
					two_xspot = None
					two_yspot = None
					two_wind = None
					two_dir = None
				else:
					# Check if two is nodata, also can't go this way
					if rastervalues[irow-1][icol+1] == -9999:
						two_xspot = None
						two_yspot = None
						two = None
						two_wind = None
						two_dir = None
					else:
						# Get wind speed resultant in two direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						two_xspot = float(xspot+cellsize)
						two_yspot = float(yspot+cellsize)
						two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]) * sqrtOf2xMidCell
						two_wind = Speed
						two_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							two_dir = cost
						else:
							two_xspot = None
							two_yspot = None
							two = None
							two_wind = None
							two_dir = None
						
			# --------------------------------
			#Position 3: Check on first column[irow][icol-1]
			#---------------------------------
			phi = 180.
			phi = phi_arr[3]
			phi_ind = 3
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
				three_wind = None
				three_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					three = None
					three_xspot = None
					three_yspot = None
					three_wind = None
					three_dir = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol-1] == -9999:
						three = None
						three_xspot = None
						three_yspot = None
						three_wind = None
						three_dir = None
					else:
						# Get wind speed resultant in zero direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						three_xspot = float(xspot-cellsize)
						three_yspot = float(yspot)
						three = (rastervalues[irow][icol]+rastervalues[irow][icol-1]) * midCell
						three_wind = Speed
						three_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							three_dir = cost
						else:
							three = None
							three_xspot = None
							three_yspot = None
							three_wind = None
							three_dir = None
						
			# -------------------------------
			#Position 4: Check on last column[irow][icol+1]
			#--------------------------------
			phi = 0.
			phi = phi_arr[4]
			phi_ind = 4
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
				four_wind = None
				four_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					four = None
					four_xspot = None
					four_yspot = None
					four_wind = None
					four_dir = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol+1] == -9999:
						four_xspot = None
						four_yspot = None
						four = None
						four_wind = None
						four_dir = None
					else:
						# Get wind speed resultant in four direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						four_xspot = float(xspot+cellsize)
						four_yspot = float(yspot)
						four = (rastervalues[irow][icol] + rastervalues[irow][icol+1]) * midCell
						four_wind = Speed
						four_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							four_dir = cost
						else:
							four_xspot = None
							four_yspot = None
							four = None
							four_wind = None
							four_dir = None
						
			# ---------------------------------------
			#Position 5: Check last row, first column[irow+1][icol-1]
			#----------------------------------------
			phi = 225
			phi = phi_arr[5]
			phi_ind = 5
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None
				five_wind = None
				five_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					five = None
					five_xspot = None
					five_yspot = None
					five_wind = None
					five_dir = None
				else:
					# Check nodata, also can't go this way
					if rastervalues[irow+1][icol-1] == -9999:
						five_xspot = None
						five_yspot = None
						five = None	
						five_wind = None
						five_dir = None
					else:
						# Get wind speed resultant in zero direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						five_xspot = float(xspot-cellsize)
						five_yspot = float(yspot-cellsize)
						five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]) * sqrtOf2xMidCell
						five_wind = Speed
						five_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							five_dir = cost
						else:
							five_xspot = None
							five_yspot = None
							five = None
							five_wind = None
							five_dir = None
						
			# ----------------------------
			#Position 6: Check on last row[irow+1][icol]
			#-----------------------------
			phi = 270.
			phi = phi_arr[6]
			phi_ind = 6
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
				six_wind = None
				six_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					six = None
					six_xspot = None
					six_yspot = None
					six_wind = None
					six_dir = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol] == -9999:
						six_xspot = None
						six_yspot = None
						six = None
						six_wind = None
						six_dir = None
					else:
						# Get wind speed resultant in six direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						six_xspot = float(xspot)
						six_yspot = float(yspot-cellsize)
						six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]) * midCell	
						six_wind = Speed
						six_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							six_dir = cost	
						else:
							six_xspot = None
							six_yspot = None
							six = None
							six_wind = None
							six_dir = None
						
			# --------------------------------------
			#Position 7: Check last column, last row[irow+1][icol+1]
			#---------------------------------------
			phi = 315.
			phi = phi_arr[7]
			phi_ind = 7
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
				seven_wind = None
				seven_dir = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					seven = None
					seven_xspot = None
					seven_yspot = None
					seven_wind = None
					seven_dir = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol+1] == -9999:
						seven = None
						seven_xspot = None
						seven_yspot = None
						seven_wind = None
					else:
						# Get wind speed resultant in seven direction
						R_x = rastervalues_u[irow][icol] + np.cos(phi)
						R_y = rastervalues_v[irow][icol] + np.sin(phi)
						Speed = np.sqrt(R_x**2 + R_y**2)
						
						seven_xspot = float(xspot+cellsize)
						seven_yspot = float(yspot-cellsize)
						seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]) * sqrtOf2xMidCell
						seven_wind = Speed
						seven_dir = wind_prob[phi_ind]
						
						cost = order_cost[np.where(phi_ind == sort_index)[0][0]]
						if not math.isnan(cost):
						#if phi_ind == closest_ind:
							seven_dir = cost
						else:
							seven = None
							seven_xspot = None
							seven_yspot = None
							seven_wind = None
							seven_dir = None
			
			# Get neighbors: x,y,wind speed			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero,zero_wind)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one,one_wind)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two,two_wind)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three,three_wind)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four,four_wind)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five,five_wind)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six,six_wind)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven,seven_wind)	
			
			'''
			# Get neighbors: x,y,cost,wind speed,direction to go			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero,zero_wind,zero_dir)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one,one_wind,one_dir)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two,two_wind,two_dir)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three,three_wind,three_dir)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four,four_wind,four_dir)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five,five_wind,five_dir)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six,six_wind,six_dir)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven,seven_wind,seven_dir)
			# Get neighbors: x,y,1-prob,wind speed,cost			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero_dir,zero_wind,zero)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one_dir,one_wind,one)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two_dir,two_wind,two)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three_dir,three_wind,three)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four_dir,four_wind,four)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five_dir,five_wind,five)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six_dir,six_wind,six)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven_dir,seven_wind,seven)
			'''
		# end:: irow loop

	# Get max and min wind speed from graph
	min_speed = np.nanmin(nWeightPairs[:,:,:,3])
	if min_speed < 0:
		pdb.set_trace() # Negative speed encountered?
	min_speed = 0
	max_speed = np.nanmax(nWeightPairs[:,:,:,3])
		
	# Convert speed to resistance
	speed_res = ((minres - maxres) / (max_speed - min_speed)) * (nWeightPairs[:,:,:,3] - max_speed) + 1
	
	# Replace wind speed with resistance values?
	nWeightPairs[:,:,:,3] = speed_res
	
	# Recalculate the cost values
	speed_cost = (nWeightPairs[:,:,:,2] / 2.) * (1 + speed_res)
	
	# Replace Euclidean cost with new wind cost
	nWeightPairs[:,:,:,2] = speed_cost
		
	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToDirection()

# --------------------------------------------------------------------------
def RasterToHikingWeights(resgrid,elevgrid,maxres,A,B,minres):
	'''
	RasterToHikingWeights()
	Create graph that is weighted based on slope, e.g., faster route around mountain
	Calculate speed
	speed = 6 * e^-3.5[slope + 0.05] is the hiking function
	speed = A * e^(B * abs([slope + C]))
	Elevation is read in and slope calculated between adjacent cells for tan(angle) = slope
	Slope values <= 0, will have value of A. 
	'''
	
	# Load resistance file
	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	
	nrows = int(header_dict['nrows'])
	ncols = int(header_dict['ncols'])
	xllcorner = float(header_dict['xllcorner'])
	yllcorner = float(header_dict['yllcorner'])
	cellsize = float(header_dict['cellsize'])
	NODATA_value = float(header_dict['NODATA_value'])
	
	# Load slope file.
	header_dict_v, data_list_v = FileIO.loadFile(elevgrid,header_lines=6)
	
	nrows_v = int(header_dict_v['nrows'])
	ncols_v = int(header_dict_v['ncols'])
	xllcorner_v = float(header_dict_v['xllcorner'])
	yllcorner_v = float(header_dict_v['yllcorner'])
	cellsize_v = float(header_dict_v['cellsize'])
	NODATA_value_v = float(header_dict_v['NODATA_value'])
	
	# Check here to make sure all the same sizes
	if nrows != nrows_v and ncols != ncols_v:
		print('Resistance and slope files are not same. (nrows x ncols)')
		sys.exit(-1)
	if xllcorner != xllcorner_v and yllcorner != yllcorner_v:
		print('Resistance and slope files are not same. (xllcorner and yllcorner)')
		sys.exit(-1)
	if cellsize != cellsize_v:
		print('Resistance and slope files are not same. (cellsize)')
		sys.exit(-1)
		
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)

	# Turn the data_list into a numpy float array
	rastervalues = np.asarray(data_list,dtype='float')
	rastervalues_v = np.asarray(data_list_v,dtype = 'float')
		
	# Leave no data as -9999, then when search around place nan there.
		
	# Sum up the resgrid cost and speed. 
	
	# Create array: x,y,resgrid cost 
	dt = np.dtype((np.float,(3,)))
	nWeightPairs = np.zeros((nrows,ncols,8),dtype=dt)
		
	# Store x and y locations of each raster pixel 
	xrasloc = []
	yrasloc = []
	
	# CSE --common-sub-expression elimination 	
	midCell = cellsize/2.0
	sqrtOf2 = np.sqrt(2.0)
	sqrtOf2xMidCell = sqrtOf2 * midCell
	# Conversion of speed to resistance
	minspeed = 0.
	maxspeed = A
	m = (maxres - minres)/(minspeed - maxspeed)
	
	# Begin irow loop
	for irow in range(nrows):
		
		yspot = yll+( cellsize * (nrows-1-irow))
		# Begin icol loop
		for icol in range(ncols):			
			
			# Get key spot name
			xspot = xll+( cellsize * icol)
			
			# Store this location value, as an integer {row,col} coordinate
			xrasloc.append(xspot)
			yrasloc.append(yspot)
				
			#Position 0: rastervalues_u[irow-1][icol-1]
			# ---------------------------------------
			# Special case for first row, on first column
			if irow == 0 or icol == 0:
				zero = None
				zero_xspot = None
				zero_yspot = None
			else:				
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					zero = None
					zero_xspot = None
					zero_yspot = None
				else:				
					# Check if zero is nodata, also can't go this way
					if rastervalues[irow-1][icol-1] == -9999:
						zero = None
						zero_xspot = None
						zero_yspot = None
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow-1][icol-1] - rastervalues_v[irow][icol]) / (sqrtOf2xMidCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						zero_xspot = float(xspot-cellsize)
						zero_yspot = float(yspot+cellsize)
						zero = (rastervalues[irow][icol]+rastervalues[irow-1][icol-1]+speed_res) * sqrtOf2xMidCell
						
			#------------------------------
			#Position 1: rastervalues_u[irow-1][icol]
			#------------------------------
			if irow == 0:
				one = None
				one_xspot = None
				one_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					one = None
					one_xspot = None
					one_yspot = None
				else:
					# Check if one is nodata, also can't go this way
					if rastervalues[irow-1][icol] == -9999:
						one = None
						one_xspot = None
						one_yspot = None
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow-1][icol] - rastervalues_v[irow][icol]) / (midCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						one_xspot = float(xspot)
						one_yspot = float(yspot+cellsize)
						one = (rastervalues[irow][icol]+rastervalues[irow-1][icol]+speed_res) * midCell
						
			# ---------------------------------------------
			#Position 2: Check on first row, on last column [irow-1][icol+1]
			#----------------------------------------------
			if irow == 0 or icol == ncols-1:
				two_xspot = None
				two_yspot = None
				two = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					two = None
					two_xspot = None
					two_yspot = None
				else:
					# Check if two is nodata, also can't go this way
					if rastervalues[irow-1][icol+1] == -9999:
						two_xspot = None
						two_yspot = None
						two = None
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow-1][icol+1] - rastervalues_v[irow][icol]) / (sqrtOf2xMidCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						two_xspot = float(xspot+cellsize)
						two_yspot = float(yspot+cellsize)
						# Add any additional speed resistance to the cost calculation
						two = (rastervalues[irow][icol]+rastervalues[irow-1][icol+1]+speed_res) * sqrtOf2xMidCell
						
			# --------------------------------
			#Position 3: Check on first column[irow][icol-1]
			#---------------------------------
			if icol == 0:
				three = None
				three_xspot = None
				three_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					three = None
					three_xspot = None
					three_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol-1] == -9999:
						three = None
						three_xspot = None
						three_yspot = None
					else:
						# Get Slope value						
						slope = math.atan((rastervalues_v[irow][icol-1] - rastervalues_v[irow][icol])/(2*midCell))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						# Add any additional speed resistance to the cost calculation
						three = (rastervalues[irow][icol]+rastervalues[irow][icol-1] + speed_res) * midCell
						three_xspot = float(xspot-cellsize)
						three_yspot = float(yspot)
						
						
			# -------------------------------
			#Position 4: Check on last column[irow][icol+1]
			#--------------------------------
			if icol == ncols-1:
				four_xspot = None
				four_yspot = None
				four = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					four = None
					four_xspot = None
					four_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow][icol+1] == -9999:
						four_xspot = None
						four_yspot = None
						four = None
					else:
						# Get Slope value						
						slope = math.atan((rastervalues_v[irow][icol+1] - rastervalues_v[irow][icol])/(2*midCell))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						# Add any additional speed resistance to the cost calculation
						four = (rastervalues[irow][icol] + rastervalues[irow][icol+1] + speed_res) * midCell
						four_xspot = float(xspot+cellsize)
						four_yspot = float(yspot)						
												
			# ---------------------------------------
			#Position 5: Check last row, first column[irow+1][icol-1]
			#----------------------------------------
			if irow == nrows-1 or icol == 0:
				five_xspot = None
				five_yspot = None
				five = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					five = None
					five_xspot = None
					five_yspot = None
				else:
					# Check nodata, also can't go this way
					if rastervalues[irow+1][icol-1] == -9999:
						five_xspot = None
						five_yspot = None
						five = None	
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow+1][icol-1] - rastervalues_v[irow][icol]) / (sqrtOf2xMidCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						five_xspot = float(xspot-cellsize)
						five_yspot = float(yspot-cellsize)
						five = (rastervalues[irow][icol] + rastervalues[irow+1][icol-1]+speed_res) * sqrtOf2xMidCell
												
			# ----------------------------
			#Position 6: Check on last row[irow+1][icol]
			#-----------------------------
			if irow == nrows-1:
				six_xspot = None
				six_yspot = None
				six = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					six = None
					six_xspot = None
					six_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol] == -9999:
						six_xspot = None
						six_yspot = None
						six = None
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow+1][icol] - rastervalues_v[irow][icol]) / (midCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						six_xspot = float(xspot)
						six_yspot = float(yspot-cellsize)
						six = (rastervalues[irow][icol] + rastervalues[irow+1][icol]+speed_res) * midCell	
						
			# --------------------------------------
			#Position 7: Check last column, last row[irow+1][icol+1]
			#---------------------------------------
			if icol == ncols-1 or irow == nrows-1:
				seven = None
				seven_xspot = None
				seven_yspot = None
			else:
				# Check if at a nodata spot, can't be here
				if rastervalues[irow][icol] == -9999:
					seven = None
					seven_xspot = None
					seven_yspot = None
				else:
					# Check if nodata, also can't go this way
					if rastervalues[irow+1][icol+1] == -9999:
						seven = None
						seven_xspot = None
						seven_yspot = None
					else:
						# Get Slope value
						slope = math.atan((rastervalues_v[irow+1][icol+1] - rastervalues_v[irow][icol]) / (sqrtOf2xMidCell * 2))
						# Get speed
						if slope < 0:
							speed = A
						else:
							speed = A * np.exp(B * slope)
						# Convert speed to resistance
						speed_res = m * speed + maxres
						seven_xspot = float(xspot+cellsize)
						seven_yspot = float(yspot-cellsize)
						seven = (rastervalues[irow][icol] + rastervalues[irow+1][icol+1]+speed_res) * sqrtOf2xMidCell
			
			# Get neighbors: x,y,wind speed			
			nWeightPairs[irow][icol][0] = (zero_xspot,zero_yspot,zero)
			nWeightPairs[irow][icol][1] = (one_xspot,one_yspot,one)
			nWeightPairs[irow][icol][2] = (two_xspot,two_yspot,two)
			nWeightPairs[irow][icol][3] = (three_xspot,three_yspot,three)
			nWeightPairs[irow][icol][4] = (four_xspot,four_yspot,four)
			nWeightPairs[irow][icol][5] = (five_xspot,five_yspot,five)
			nWeightPairs[irow][icol][6] = (six_xspot,six_yspot,six)
			nWeightPairs[irow][icol][7] = (seven_xspot,seven_yspot,seven)	
		
		# end:: irow loop
	
	return nWeightPairs, xrasloc,yrasloc, header_dict
	
	# End::RasterToHikingWeights()
	
	
if scipyAvail:	
	# --------------------------------------------------------------------------
	def MatchToRaster(xyfilename,xrasloc,yrasloc,directionans):
		'''
		MatchToRaster() Pre XY manipulation, read in and assign key value for NWeights 
		'''
		header_dict, data_list = FileIO.loadFile(xyfilename,header_lines=1,delimiter=',')
		xy_points = np.asarray(data_list,dtype='float')
		# Extract values, make numpy
		data = np.zeros((len(xrasloc),2))
		data[:,0] = xrasloc
		data[:,1] = yrasloc
		
		# Call the KDTree function to get closest values
		tree = KDTree(data)
		
		# Query the tree to get fixed points
		fixed_pts = tree.query(xy_points)
		
		# The re-write to match G convention
		stringpts = []
		for i in range(len(fixed_pts[1])):
			#s1 = str(float(xrasloc[fixed_pts[1][i]]))
			#s2 = str(float(yrasloc[fixed_pts[1][i]]))
			s1 = str(tree.data[fixed_pts[1][i]][0])
			s2 = str(tree.data[fixed_pts[1][i]][1])
			
			if len(stringpts) != 0:
			
				if (np.asarray(stringpts) == s1 + '_' + s2).any():
					
					print(('Warning: There are overlapping points around x = ' + s1 + ' y = ' + s2 + '.'))
					if directionans:
						print('Overlapping points must be removed for directional answer TRUE.')
						sys.exit(-1)					
					
				stringpts.append(s1+ '_' + s2)
				
			else:
				# concatenate two additive terms...
				stringpts.append(s1+ '_' + s2)
				# end:: i fixed_pts[] loop
		
		# Return values
		tupMtoG = stringpts, len(stringpts)
		return tupMtoG
			
		# End::MatchToRaster()

# ------------------------------------------------------------------------------------------
def isApprox(a,b,TOLERANCE = 10E-6):
	'''
	isApprox() performs a boolean comparison of 2 floating 
	point variables in a safe manner. For a to 'equal' b, they simply 
	must be within TOLERANCE distance of each other. Be careful to set
	TOLERANCE with awareness of single vs. double precision of the
	quantities compared.
	'''
	rsltBool = False
	if math.fabs(a-b) <= TOLERANCE:
		rsltBool = True
	# return True if a is within TOLERANCE of b, otherwise False
	return rsltBool
		
# ------------------------------------------------------------------------------------------
def calcOutputs(indexDict,isource,idest,path,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathlen):
	'''
	calcOutputs()
	Calculates path additions and cost distance matrix
	'''
	
	while True:
	
		# If no path calculated
		if len(path) == 0:
			
			for i in range(len(idest)):	
				# Then write shortest path len value to matrix
				cdmatrix[indexDict[isource]][indexDict[idest[i]]] = None
				
				# Append path information
				paths.append((isource,idest[i],None,None))
		
		# If None path calculated
		elif len(path) != 0 and path[0] == None:
			
			# Then write shortest path len value to matrix
			cdmatrix[indexDict[isource]][indexDict[idest]] = None
				
			# Append path information
			paths.append((isource,idest,None,None))
				
		# If path calculated
		elif len(path) != 0 and path[0] != None:
			
			for ipath in range(len(path)):
							
				(xspot,yspot) = setFloat(path[ipath])
				
				x = approxWhere(xvalues,xspot,10E-3)
				y = approxWhere(yvalues,yspot,10E-3)				
							
				# sum to the x,y coordinate of the new path piece
				pathadd[y[0][0]][x[0][0]] += 1				
						
			# Append path information
			paths.append((isource,idest,pathlen,path))        
			
			# Then write shortest path len value to matrix
			cdmatrix[indexDict[isource]][indexDict[idest]] = pathlen
			
		# return desired output data
		return pathadd, cdmatrix, paths

		# End::calcOutputs()

		# ------------------------------------------------------------------------------------------
def calcOutputs_wind(indexDict,isource,idest,path,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathlen,pathcond,cdmatrix_accCost):
	'''
	calcOutputs()
	Calculates path additions and cost distance matrix
	pathlen is wind direction
	pathcond is cost LCP
	'''
	
	while True:
	
		# If no path calculated
		if len(path) == 0:
			
			for i in range(len(idest)):	
				# Then write shortest path len value to matrix
				cdmatrix[indexDict[isource]][indexDict[idest[i]]] = None
				cdmatrix_accCost[indexDict[isource]][indexDict[idest[i]]] =None
				
				# Append path information
				paths.append((isource,idest[i],None,None,None))
		
		# If None path calculated
		elif len(path) != 0 and path[0] == None:
			
			# Then write shortest path len value to matrix
			cdmatrix[indexDict[isource]][indexDict[idest]] = None
			cdmatrix_accCost[indexDict[isource]][indexDict[idest]] =None
				
			# Append path information
			paths.append((isource,idest,None,None,None))
				
		# If path calculated
		elif len(path) != 0 and path[0] != None:
			
			for ipath in range(len(path)):
							
				(xspot,yspot) = setFloat(path[ipath])
				
				x = approxWhere(xvalues,xspot,10E-3)
				y = approxWhere(yvalues,yspot,10E-3)				
							
				# sum to the x,y coordinate of the new path piece
				pathadd[y[0][0]][x[0][0]] += 1				
						
			# Append path information
			paths.append((isource,idest,pathlen,pathcond,path)) 
			
			# Then write shortest path len value to matrix
			cdmatrix[indexDict[isource]][indexDict[idest]] = pathlen
			cdmatrix_accCost[indexDict[isource]][indexDict[idest]] = pathcond
			
			
		# return desired output data
		return pathadd, cdmatrix, paths,cdmatrix_accCost

		# End::calcOutputs_wind()

		
# --------------------------------------------------------------------------
def do_work(q,out,nWeightPairs,edge_type,edge_dist,xvalues,\
			yvalues):
	'''
	do_work()
	Call shortestPath functions which calls
	Dijkstra's Algorithm.
	'''         
	while True:
		#info('do_work')
		sourcept,destpts = q.get()
		if edge_type == 'normal':
			paths,vOrder,pathlens = Dijkstra(nWeightPairs,xvalues,\
			yvalues,sourcept,destpts)
		elif edge_type == 'threshold':
			paths,vOrder,pathlens =Dijkstrathreshold(edge_dist,nWeightPairs,xvalues,\
			yvalues,sourcept,destpts)
		elif edge_type == 'all_paths':
			all_paths_lens =DijkstraAllPaths(edge_dist,nWeightPairs,xvalues,\
			yvalues,sourcept)
	
		if edge_type == 'all_paths':
			out.put((all_paths_lens))
		else:
			i = 0
			
			for p in paths:
				if (p[0] != None):
					out.put((sourcept,vOrder[i],p,pathlens[i]))
				i += 1
		q.task_done()
		# End::do_work()

#---------------------------------------------------------------------------
def calc_points(nopoints,stringpts,visited,EDthresholdans=False,\
nbhd_dist=None):
	'''
	calc_points() - create the source point list (srclist) and the
	destination point list (destlist) to be used to run
	the shortest path algorithm
	'''
	
	destlist = []
	srclist = []
	indexDict = {}
	# Begin source loop
	for isource in range(nopoints):
		idestpts = []
	 	
		indexDict[stringpts[isource]] = isource
		#add the source point to a dict and associated index value to later create the 
		# cdmatrix
		
		# Start destination loop
		for idest in range(nopoints):
		  	#add the dest point to a dict and associated index value to later create the 
			# cdmatrix
			# Check for looping backwards and if itself
			if isource == idest or visited[isource][idest] == 1:
				break
			
			# The case where ED between points is a threshold
			elif EDthresholdans:
				
				isource_x = float(stringpts[isource].split('_')[0])
				isource_y = float(stringpts[isource].split('_')[1])
				idest_x = float(stringpts[idest].split('_')[0])
				idest_y = float(stringpts[idest].split('_')[1])
	
				sourcedist_to_dest = np.sqrt((isource_x-idest_x) **2.0+\
				(isource_y-idest_y) **2.0)
				
				# Then add the paired_pts to calculate
				if sourcedist_to_dest < nbhd_dist:
					
					idestpts.append(stringpts[idest])
			# Else grab all points       
			else:
				
				idestpts.append(stringpts[idest])
			# Store visited paths
			visited[isource][idest] = 1
		if len(idestpts) > 0:
			destlist.append(idestpts)

			srclist.append(stringpts[isource])
	
	# Return value
	return srclist, destlist, indexDict

#---------------------------------------------------------------------------
def accumulateDistances(trans_func,full_paths,pathadd,xvalues,yvalues,srclist,const_kernal_vol,edge_dist,vol_constant): 
	count = 0
	for dist_dicts in full_paths:
		temp_distances = np.zeros((pathadd.shape))
		for point,distance in dist_dicts.items():
			
			(xspot,yspot) = setFloat(point)
		
			x = approxWhere(xvalues,xspot,10E-3)
			y = approxWhere(yvalues,yspot,10E-3)				
			#if x == 53 and y == 44:
			# sum to the x,y coordinate of the new path piece
			temp_distances[y[0][0]][x[0][0]] = distance			
				
		if trans_func == 'linear':
			temp_distances = 1. - (1./edge_dist) * temp_distances
			# Set values outside of threshold back to 0
			temp_distances[temp_distances == 1.0] = 0.0
			# But not the source point
			# Add code here to fix that point
			
		
		#temp_distances = temp_distances/np.max(temp_distances)
		if trans_func == 'inverse_square':
			temp_distances = 1./(temp_distances**2. + 1)
			temp_distances[temp_distances == np.inf] = np.nan
			# Set values outside of threshold back to 0
			temp_distances[temp_distances == 1.0] = 0.0
			# But not the source point
			# add code here to fix this
		
		if const_kernal_vol:		
			vol_const = vol_constant * 3/(math.pi*edge_dist**2)
		else:
			vol_const = 1.
		
		pathadd += temp_distances * vol_const 
	

	'''
	for pts in srclist:
		isource_x = float(pts.split('_')[0])
		isource_y = float(pts.split('_')[1])
		
		x = approxWhere(xvalues,isource_x,10E-3)
		y = approxWhere(yvalues,isource_y,10E-3)				
		pathadd[y[0][0]][x[0][0]] = len(srclist) + 1	
	'''
	return pathadd
		
#---------------------------------------------------------------------------
def parallel_paths(logfHndl,nopoints,stringpts,visited,\
	EDthresholdans,num_of_pro,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,\
	edge_type,edge_dist,trans_func,const_kernal_vol,vol_constant = None,nbhd_dist=None,):  
	
	'''
	parallel_paths()
	If Linux operating system is used, then run user specified number of
	processes.
	'''
	# The queue containing the sets of source and destination points 
	# for shortest path function        
	work_queue = JoinableQueue()
	
	# Contains the output paths as individual lists of single paths
	out_queue = JoinableQueue()
	
	start_time = datetime.datetime.now()
	
	# Create the process list
	processes = [Process(target=do_work, \
	args=(work_queue,out_queue,nWeightPairs,edge_type,edge_dist,xvalues,\
			yvalues)) for i in range(num_of_pro)]    

	# Start the processes
	for p in processes:
		p.daemon = True
		p.start()
	
	# Calc the points to be used in the shortest path algorithm
	srclist, destlist, indexDict  = calc_points(nopoints,stringpts,visited,\
	EDthresholdans,nbhd_dist)    
	
	if edge_type == 'all_paths':
		for i, srcpt in enumerate(stringpts):
			work_queue.put((srcpt,0.))
	#loop through the source and destinations and put them on the work_queue
	else:
		for i in range(len(destlist)):
			 work_queue.put((srclist[i], destlist[i]))
	
	# Blocks until all processes are completed
	work_queue.join()       
	
	#Create a list to append values from the out_queue to then calculate outputs
	pathList = []
	full_paths = []
	

	stringout = '\nTotal time path calculation: '+\
	str(datetime.datetime.now() - start_time)+'\n'

	FileIO.logMsg(logfHndl,stringout)
	start_time = datetime.datetime.now()
	
	# Loop to do the output calculations
	for i in range(out_queue.qsize()):
		
		if edge_type == 'all_paths':
			all_paths_len = out_queue.get()
			full_paths.append(all_paths_len)
		else:				
			isource,idest,path,pathlen = out_queue.get()
			pathList.append((isource,idest,path,pathlen))	
		
	if edge_type == 'all_paths':
		pathadd = accumulateDistances(trans_func,full_paths,pathadd,xvalues,yvalues,srclist,const_kernal_vol,edge_dist,vol_constant) 
		tupCalcOut = (pathadd,0,0)
	
	else:
		if len(pathList) > 0:
			for i in range(len(pathList)):
				isource,idest,path,pathlen = pathList[i]
				
				tupCalcOut = calcOutputs(indexDict,isource,idest,path,nWeightPairs,\
				xvalues,yvalues,pathadd,cdmatrix,paths,pathlen)
		else:
			print("The path list is empty, try increasing the path threshold")
	
	stringout = '\nTotal time for calculation outputs:  '+\
	str(datetime.datetime.now() - start_time)+'\n'
	FileIO.logMsg(logfHndl,stringout)
	#tupCalcOut = (pathadd,cdmatrix,paths)				
	for p in processes:
		p.terminate()	
	return tupCalcOut 

	# End::parallel_paths()

#---------------------------------------------------------------------------
def serial_paths(dirtype,logfHndl,nopoints,stringpts,visited,EDthresholdans,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,edge_type,edge_dist,trans_func,const_kernal_vol,vol_constant =None,nbhd_dist=None,resans=None,cdmatrix_accCost=None):
	'''
	serial_paths() If windows operating system then run this function.
	'''

	# Out Queue
	out_queue = queue.Queue(0)
	
	# Call function to calculate the points - lower
	srclist, destlist, indexDict = calc_points(nopoints,stringpts,visited,\
	EDthresholdans,nbhd_dist)  
	
	# Timing events: start and log print out
	start_time = datetime.datetime.now()	
	spathsprint = 'Printing out every source,destination combination at 10 percent completetions...'
	FileIO.logMsg(logfHndl,spathsprint)
	# Count variable for log file
	count = 1	
	
	# For the length of paired points
	for i in range(len(srclist)): 
		
		# Print to log file only every 10th
		if np.mod(i,round(len(srclist)/10.)) == 0.0:
			spathsprint = "Percent complete: "+str(count*10)+"%; Source point: "+str(srclist[i])+"; Start time: "+\
			str(datetime.datetime.now()-start_time)
			FileIO.logMsg(logfHndl,spathsprint)
			count +=1
		
		# If threshold specified as 'max'
		if edge_type == 'normal':
			if resans: # Use resistance LCP
				#if dirtype != 'Wind':
				path, vOrder,pathlens =Dijkstra(nWeightPairs,xvalues,yvalues,srclist[i],destlist[i])
				#else: # pathlens is wind cost, pathcond is the cost along wind path
				#	path, vOrder,pathlens,pathcond =Dijkstra_wind(nWeightPairs,xvalues,\
				#	yvalues,srclist[i],destlist[i])
			else: # Use conductance
				
				path, vOrder,pathlens,pathcond =Dijkstra_cond(nWeightPairs,xvalues,\
				yvalues,srclist[i],destlist[i])
		# If threshold specified as not 'max'
		elif edge_type == 'threshold':
			path, vOrder,pathlens =Dijkstrathreshold(edge_dist,nWeightPairs,xvalues,\
			yvalues,srclist[i],destlist[i])	
		elif edge_type == 'all_paths':
			path =DijkstraAllPaths(edge_dist,nWeightPairs,xvalues,\
			yvalues,srclist[i])
		
		if (len(path) > 0):	
			if edge_type == 'all_paths':
				# need to pack the path into a list for accumulateDistance function which is
				# expecting a list of Dict items
				pathadd = accumulateDistances(trans_func,[path],pathadd,xvalues,yvalues,srclist,const_kernal_vol,edge_dist,vol_constant) 
				tupCalcOut = (pathadd,0,0)
			else:
				
				for j in range(len(vOrder)):
					
					if resans: # Use resistance
						tupCalcOut = calcOutputs(indexDict,srclist[i],vOrder[j],path[j],nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathlens[j])
						
					else: # Use conductance
						tupCalcOut = calcOutputs(indexDict,srclist[i],vOrder[j],path[j],nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathcond[j])
		else: 
			
			print("The path list is empty, try increasing the path threshold, caclulating path as Nan")
			
			
			if resans: # Use resistance
				
				#if dirtype == 'Wind':
				#tupCalcOut = calcOutputs_wind(indexDict,srclist[i],destlist[i],path,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathlens,pathcond,cdmatix_accCost)
				#else:
				tupCalcOut = calcOutputs(indexDict,srclist[i],destlist[i],path,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathlens)
			else: # Use conductance
				tupCalcOut = calcOutputs(indexDict,srclist[i],destlist[i],path,nWeightPairs,xvalues,yvalues,pathadd,cdmatrix,paths,pathcond)
	
	# Return values
	return tupCalcOut
	
	# End::serial_paths
	