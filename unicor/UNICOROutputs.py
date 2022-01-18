# ---------------------------------------------------------------------------------------------------
# UNICOROutputs.py
# Author(s): Erin L Landguth, Brian Hand, Joe Glassy
# revised: 2022-01-18 ELL 
# Created: September 2009
# v 0.5 -- 2010 August 30: Code Cleanup and function additions to outputs.
# ------------------------------------------------------------------------------------------------------
# Messy, fix later
msgVerbose = False
#outputGraphMetrics = False

# Import Modules with Except/Try statements

# For ELL debugging
import pdb

# Numpy functions
try:
    import numpy as np                    
except ImportError:
    raise ImportError("Numpy required.")

# Timing and general functions
try:
	import time, datetime, copy, sys, os
except ImportError:
	raise ImportError("Time and Datetime and copy required.")

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

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ------------------------------------------------------------------------------------------
def output_cdmatrix(logfHndl,cd_matrix, output_filename,data_dir,stringpts): 
	'''
	output_cdmatrix()
	This writes the cost distance matrix to file.
	Checks for duplicate points
	'''
	
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Create a full matrix, recall that cdmatrix comes in as a lower tri-diagonal
	temp = np.transpose(cd_matrix)
	cd_matrix = temp + cd_matrix
	
	# Check for duplicates and copy lines
	uni_strpts = count_unique(stringpts)
	duplicates_strpts = np.where(uni_strpts[1]>1)[0]
	if len(duplicates_strpts) != 0:
		# Loop through duplicates
		for i in range(len(duplicates_strpts)):
			
			# Get dupliccate value
			strdouble = uni_strpts[0][duplicates_strpts[i]]
			# Find location in strpts
			strpts_loc = np.where(strdouble == np.asarray(stringpts))[0]
			# Then copy the last duplicate (always the last) into the other lines
			copyrow = cd_matrix[strpts_loc[-1],:]
			copycol = cd_matrix[:,strpts_loc[-1]]
			for jspot in strpts_loc:
				cd_matrix[jspot,:] = copyrow
				cd_matrix[:,jspot] = copycol
				# Then make sure that all dubplicate spots to cuplicate spots are 0
				for kspot in strpts_loc:
					cd_matrix[jspot][kspot] = 0.
					
	# Make sure diagnols are 0
	cd_matrix[np.diag_indices_from(cd_matrix)] = 0.		
		
	FileIO.outputGrid(data_dir+output_filename,cd_matrix,delimiter=',')
	
	stringout = 'The file '+str(data_dir+output_filename)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')

	# End::output_cdmatrix()

# ------------------------------------------------------------------------------------------
def output_cdmatrix_asymmetric(logfHndl,cd_matrix_lower,cd_matrix_upper, output_filename,data_dir,stringpts,resans): 
	'''
	output_cdmatrix_asymetric()
	This writes the cost distance matrix to file for asymmetric costs.
	lower is ordered 1,2. upper is reverse 2,1.
	Use lower's lower triangle and fill lower cdmatrix.
	Use upper's lower triangle and fill upper cdmatrix.
	
	'''
	
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Create a full matrix, recall that cdmatrix comes in as a lower tri-diagonal
		
	# Matrices are read by row, want to flip them to read by column: upper is high to low ordering so different flipping needed.
	
	# Lower transpose
	cd_matrix_lower = np.transpose(cd_matrix_lower) # This makes it upper
	# Then make lowertri here all 0
	cd_matrix_lower[np.tril_indices_from(cd_matrix_lower)] = 0.
	
	# Upper want lr, tran, lr
	cd_matrix_upper = np.fliplr(np.transpose(np.fliplr(cd_matrix_upper))) # This makes it lower
	# Then make uppertri here all 0
	cd_matrix_upper[np.triu_indices_from(cd_matrix_upper)] = 0.
	
	# Then add
	cd_matrix = cd_matrix_lower + cd_matrix_upper 
	
	# If res, make diagnols 0
	if resans:
		cd_matrix[np.diag_indices_from(cd_matrix)] = 0.
	# Else, make 1
	else:	
		cd_matrix[np.diag_indices_from(cd_matrix)] = 1.
	
	'''
	# Lower ordered low to high, so different flip needed
	cd_matrix_lower = np.flipud(np.transpose(np.fliplr(cd_matrix_lower)))
	# Set diagnols and upper lower to zeros
	cd_matrix_lower[np.diag_indices_from(cd_matrix_lower)] = 0.
	cd_matrix_lower[np.tril_indices_from(cd_matrix_lower)] = 0.
	
	# Upper flip and use for the lower triangle (i know confusing)
	cd_matrix_upper = np.fliplr(np.transpose(np.fliplr(cd_matrix_upper)))
	# Set diagnols and upper triangle to zeros
	#if resans:
	cd_matrix_upper[np.diag_indices_from(cd_matrix_upper)] = 0.
	cd_matrix_upper[np.triu_indices_from(cd_matrix_upper)] = 0.
	#else:
		#cd_matrix_upper[np.diag_indices_from(cd_matrix_upper)] = 1.
	#	cd_matrix_upper[np.triu_indices_from(cd_matrix_upper)] = 0.
	
	cd_matrix = cd_matrix_lower + cd_matrix_upper 
	if resans: # Resistance
		cd_matrix[np.diag_indices_from(cd_matrix)] = 0.
	else: # Conductance		
		# and convert nans to 0
		cd_matrix = np.array(cd_matrix,dtype=str)
		cd_matrix[np.where(cd_matrix == 'nan')[0]] = '0.0'
		cd_matrix = np.array(cd_matrix,dtype=float)
		cd_matrix[np.diag_indices_from(cd_matrix)] = 1.
	'''
	# Write to file
	FileIO.outputGrid(data_dir+output_filename,cd_matrix,delimiter=',')
	
	stringout = 'The file '+str(data_dir+output_filename)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')

	# End::output_cdmatrix_asymmetrical()
	
	
# ------------------------------------------------------------------------------------------
def output_pathadd(logfHndl,path_values,pathadd_file_name,header_dict,keys,data_dir):
	'''
	output_pathadd()
	Writes path additions to ascii file with header information
	'''    
	
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Turn 0.0 pathadd values back to nodata -9999 values
	#path_values[path_values==0.0] = -9999.0
	
	FileIO.outputGrid(data_dir+pathadd_file_name,path_values,header_dict = header_dict,keys = keys)
	
	stringout = 'The file '+str(data_dir+pathadd_file_name)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')

	# End::output_pathadd()
	
#-----------------------------------------------------------------------------------
def output_paths(logfHndl,paths,paths_file_name,data_dir):
	'''
	output_paths() - output the full list of paths between points, 
	and shortest path distance for a given path
	'''
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Open file to write to
	fout = open(data_dir+paths_file_name, 'w')
	
	# Write out title information
	fout.write('Source,Destination,PathLength,PathConnections\n')
	
	# Begin loop through point combinations
	for ipath in range(len(paths)):
	
		# Print out path source, destination, length
		fout.write(paths[ipath][0])
		fout.write(',')
		fout.write(paths[ipath][1])
		fout.write(',')
		fout.write(str(paths[ipath][2]))
		fout.write(',')
		# Print out the path connections
		if paths[ipath][3] != None:
			for ipathlen in range(len(paths[ipath][3])):
				fout.write(paths[ipath][3][ipathlen])
				fout.write(',')
		else:
			fout.write('nan')
		fout.write('\n')
		
	# Logging information
	stringout = 'The file '+str(data_dir+paths_file_name)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')
	
	# Close file
	fout.close
	# End::output_paths
	
# ------------------------------------------------------------------------------------------
def output_buffer(logfHndl,buff_values,buff_file_name,header_dict,keys,data_dir):
	'''
	output_buffer()
	Writes the buffered path values to ascii file with header information
	'''  

	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Turn 0.0 pathadd values back to nodata -9999 values
	buff_values[buff_values==0.0] = -9999.0
	
	FileIO.outputGrid(data_dir+buff_file_name,buff_values,header_dict = header_dict,keys = keys)
	
	stringout = 'The file '+str(data_dir+buff_file_name)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')

	# End::output_buffer()
	
def output_levels(logfHndl,level_values,level_file_name,header_dict,keys,data_dir):
	'''
	output_levels()
	Writes the catagorical buffered path values to ascii file with header information
	'''    
	
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	FileIO.outputGrid(data_dir+level_file_name,level_values,header_dict = header_dict,keys = keys)

	stringout = 'The file '+str(data_dir+level_file_name)+' has been created in: '+\
	str(datetime.datetime.now() -start_time1)
	FileIO.logMsg(logfHndl,stringout+'\n')	

	# End::output_levles()		

#if networkxAvail:	
	# ------------------------------------------------------------------------------------------
def output_graphmetrics(pathadd,paths,file_name,data_dir):
	'''
	output_graphmetrics()
	Calculates graph theory metrics from package NetworkX, stores in file and
	outputs in .csv file.
	''' 
	
	# Graph package
#print outputGraphMetrics
#if outputGraphMetrics:
	try:
		import networkx	as nx	
		
	except ImportError:
		
		raise ImportError("NetworkX required.")

	pathG = nx.Graph()
	
	# Get nrows and columns
	nrows = len(pathadd)
	ncols = len(pathadd[0])
	
	# Now loop through pathadd rows 
	for irow in range(nrows):
		
		# Begin loop through pathadd columns
		for icol in range(ncols):
			
			# Skip -9999. values
			if pathadd[irow][icol] != -9999.0:
				
				# Get spot node number 
				nodenumber = ncols*irow+icol
								
				# Add node to pathG
				pathG.add_node(nodenumber)
				
				# Check neighbors for edges all 8 conditions
				# Count 0,1,2,3,4,5,6,7 in counter-clockwise
				# around center cell starting 0 in lower
				# left corner
				
				# Left top corner:
				if irow == 0 and icol == 0:
				
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 2
					if pathadd[irow+1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
				
				# Right top corner
				elif irow == 0 and icol == ncols-1:
					
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 0
					if pathadd[irow+1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
				
				# Left bottom corner
				elif irow == nrows-1 and icol == 0:
				
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 4
					if pathadd[irow-1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
				
				# Right bottom corner
				elif irow == nrows-1 and icol == ncols-1:
					
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
						
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
						
					# Get neighbors: spot 6
					if pathadd[irow-1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
						
				# Top side
				elif irow == 0 and icol != 0 and icol != ncols-1:
				
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 0
					if pathadd[irow+1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 2
					if pathadd[irow+1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
				
				# Left side
				elif icol == 0 and irow != 0 and irow != nrows-1:
				
					# Get neighbors -spots 1,2,3,4,5
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 2
					if pathadd[irow+1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 4
					if pathadd[irow-1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
									
				# Right side
				elif icol == ncols-1 and irow != 0 and irow != nrows-1:
				
					# Get neighbors - spots 0,1,5,6,7
					# Get neighbors: spot 0
					if pathadd[irow+1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 6
					if pathadd[irow-1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
										
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
												
				# Bottom side:
				elif irow == nrows-1 and icol != 0 and icol != ncols-1:
				
					# Get neighbors - spots 3,4,5,6,7
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 4
					if pathadd[irow-1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
											
					# Get neighbors: spot 6
					if pathadd[irow-1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
										
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
													
				# Everything else:
				else:
				
					# Get neighbors: spot 0
					if pathadd[irow+1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 1
					if pathadd[irow+1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 2
					if pathadd[irow+1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow+1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 3
					if pathadd[irow][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 4
					if pathadd[irow-1][icol+1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol+1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
					
					# Get neighbors: spot 5
					if pathadd[irow-1][icol] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+icol
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
										
					# Get neighbors: spot 6
					if pathadd[irow-1][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow-1)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
										
					# Get neighbors: spot 7
					if pathadd[irow][icol-1] != -9999.:
						
						# Then get egde number
						edgenumber = ncols*(irow)+(icol-1)
						
						# Then add edge to pathG
						pathG.add_edge(nodenumber,edgenumber)
	
	# Calculate properties from path lengths: min, max, average
	pathlen = []
	for i in range(len(paths)):
		pathlen.append(paths[i][2])
	
	# Create file to write info to
	try:
		fout = open(data_dir+file_name, 'w')
	except(IOerror,OSerror) as e:
		print("UNICOROutputs  %s, error%s",(filename,e))
		sys.exit(-1)
		
	# Write header information
	fout.write('Minimum Path Length,')
	fout.write(str(min(pathlen))+'\n')
	fout.write('Maximum Path Length,')
	fout.write(str(max(pathlen))+'\n')
	fout.write('Average Path Length,')
	fout.write(str(sum(pathlen)/len(paths))+'\n')
	fout.write('Density of Graph,')
	fout.write(str(nx.density(pathG))+'\n')
	fout.write('Number of nodes,')
	fout.write(str(nx.number_of_nodes(pathG))+'\n')
	fout.write('Number of edges,')
	fout.write(str(nx.number_of_edges(pathG))+'\n')
	fout.write('Is the graph a bipartite,')
	fout.write(str(nx.is_bipartite(pathG))+'\n')
	fout.write('Size of the largest clique,')
	fout.write(str(nx.graph_clique_number(pathG))+'\n')
	fout.write('Number of maximal cliques,')
	fout.write(str(nx.graph_number_of_cliques(pathG))+'\n')
	fout.write('Transitivity,')
	fout.write(str(nx.transitivity(pathG))+'\n')
	fout.write('Average clustering coefficient,')
	fout.write(str(nx.average_clustering(pathG))+'\n')
	fout.write('Test graph connectivity,')
	fout.write(str(nx.is_connected(pathG))+'\n')
	fout.write('Number of connected components,')
	fout.write(str(nx.number_connected_components(pathG))+'\n')
	fout.write('Consists of a single attracting component,')
	fout.write(str(nx.is_attracting_component(pathG))+'\n')
	if nx.is_attracting_component(pathG) == True:
		fout.write('Number of attracting components,')
		fout.write(str(nx.number_attracting_components(pathG))+'\n')
	if nx.is_connected(pathG):
		fout.write('Center,')
		fout.write(str(nx.center(pathG))+'\n')
		fout.write('Diameter,')
		fout.write(str(nx.diameter(pathG))+'\n')
		#fout.write('Eccentricity,')
		#fout.write(str(nx.eccentricity(pathG))+'\n')
		fout.write('Periphery,')
		fout.write(str(nx.periphery(pathG))+'\n')
		fout.write('Radius,')
		fout.write(str(nx.radius(pathG))+'\n')
	fout.write('Degree assortativity,')
	fout.write(str(nx.degree_assortativity(pathG))+'\n')
	fout.write('Degree assortativity Pearsons r,')
	fout.write(str(nx.degree_pearsonr(pathG))+'\n')
			
	# Close file
	fout.close
	
	del(pathG)
	# End::output_graphmetrics()
		