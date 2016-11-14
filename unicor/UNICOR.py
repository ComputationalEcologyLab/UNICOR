# -------------------------------------------------------------------------------------------------------------
# UNICOR.py
# Author(s): Erin L Landguth, Brian Hand, Joe Glassy
# most recent revision: "2011-03-07T12:07:00MST"(jmg)
# originally created  : September 2009
# ------------------------------------------------------------------------------------------------------
appRele = "2015-10-20T13:20:00MDT"
appName = "UNICOR"
# jmg corrected appVers tag to reflect actual software version (not SVNs, and not pending '1.0')
appVers = "2.5"
authorNames = "E.L. Landguth, B.K. Hand, J.M. Glassy"

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# Point threshold for graph metrics NetworkX library
noptsthreshold = 100
DEBUG_MODE = False
#File absolute paths for importing functions
UTILITIES_PATH =  "../utilities/"
keys = ['ncols','nrows','xllcorner','yllcorner',\
		'cellsize','NODATA_value']

# Import Modules with Except/Try statements

# For ELL debugging
import pdb

# Platform and system functions
try:
	import os, sys                    
except ImportError as eMsg:
	print("ImportError: (%s) OS and SYS required!"%(eMsg))
	import sys
	sys.exit(-1)

#Import the package specific folders
RIPMGR_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"RipMgr"))
FILEIO_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"FileIO"))
DEBUG_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"DebugTools"))

if RIPMGR_folder not in sys.path:
     sys.path.insert(0, RIPMGR_folder)
if FILEIO_folder not in sys.path:
     sys.path.insert(0, FILEIO_folder)
if DEBUG_folder not in sys.path:
     sys.path.insert(0, DEBUG_folder)
	
if DEBUG_MODE:
	import DebugTools as DBug

# Numpy functions
try:
	import numpy as np                    
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)		

# UNICOR functions
try:
	from UNICORFun import *    
	from UNICORMaps import *
	from UNICOROutputs import *
except ImportError as eMsg:
	print("ImportError: (%s) UNICORFun, UNICORMaps, and/or UNICOROutputs required."%(eMsg))
	print("NOTE: failure of deeper dependencies (missing scipy,networkx etc) can also trigger this error!")
	sys.exit(-1)
	
# Utility Functions
try:
	from RipMgr import mgrParseStanzaInputs    # Joe's Rip Manager for inputs
	from RipMgr import *
	import FileIO
except ImportError as eMsg:
	print("ImportError: RipMgr module(from Lupine Logic Inc) and FileIO.py are required."%(eMsg)) 

def main(ripFilePath):
#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 

	# Timing events: start
	start_time0 = datetime.datetime.now()		
	
	# If .rip file does not exist
	if not os.path.exists(ripFilePath):
		print("Cannot find or open runtime inputs file(%s)"%(ripFilePath))
		sys.exit(-1)
	
	# --------------------------------------------------
	# Get user defined program input
	# --------------------------------------------------
			
	# create a RipMgr instance, via parser itself
	r = mgrParseStanzaInputs(ripFilePath)   
	sessionLbl = r.kwdGetValue('Session_label')
	# This properly names log file
	sessionLbl = sessionLbl.strip(' ')
	logSessionPath = sessionLbl + ".log"
	try:
		logfHndl =open(logSessionPath,'w')
	except (IOError,OSError) as eMsg:
		print("Error (%s) opening session logfile(%s)"%(eMsg,logSessionPath))
		sys.exit(-1)
	
	msgVerbose = True	
	FileIO.logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers),msgVerbose)
	FileIO.logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n',msgVerbose)
	FileIO.logMsg(logfHndl,"Session runtime inputs from: (%s)"%(ripFilePath)+'\n\n',msgVerbose)   
	FileIO.logMsg(logfHndl,"Log output directed to     : (%s)"%(logSessionPath),msgVerbose)
	#---------------------------------------------------------------------------
	# Start main     
	#---------------------------------------------------------------------------
	
	gridfile = r.kwdGetValue('Grid_Filename')
	xyfilename = r.kwdGetValue('XY_Filename')
	resans = r.kwdGetValue('Use_Resistance')
	directionans = r.kwdGetValue('Use_Direction')
	dirtype = r.kwdGetValue('Type_Direction')
	barrfile = r.kwdGetValue('Barrier_or_U_Filename')
	elevfile = r.kwdGetValue('Direction_or_V_Filename')
	minmaxres = r.kwdGetValue('Speed_To_Resistance_Scale')
	EDthresholdans = r.kwdGetValue('Use_ED_threshold')
	nbhd_dist      = float(r.kwdGetValue('ED_Distance'))
	edge_dist      = float(r.kwdGetValue('Edge_Distance'))
	edge_type = r.kwdGetValue('Edge_Type')
	num_of_pro     = int(r.kwdGetValue('Number_of_Processes'))
	outputPathadd  = r.kwdGetValue('Save_Path_Output')
	outputPaths    = r.kwdGetValue('Save_IndividualPaths_Output')
	outputGraphMetrics = r.kwdGetValue('Save_GraphMetrics_Output')
	KernelFunction = r.kwdGetValue('KDE_Function')
	BufferGridSize = r.kwdGetValue('KDE_GridSize')
	outputBuffer   = r.kwdGetValue('Save_KDE_Output')
	LevelNumber    = int(r.kwdGetValue('Number_of_Categories'))
	outputLevels   = r.kwdGetValue('Save_Category_Output')
	CDmatrixans    = r.kwdGetValue('Save_CDmatrix_Output') 
	transform_func = r.kwdGetValue('Transform_function')	
	const_kernal_vol = r.kwdGetValue('Const_kernal_vol')
	vol_constant = int(r.kwdGetValue('Kernel_volume'))

	# -----------------------------------------
	# Some error checking with mutliple options
	# -----------------------------------------
	# For the different directional models
	if directionans:
		tempdirtype = dirtype.split(';')
		if len(tempdirtype) > 1: # This is the hiking application
			if not resans:
				print('Resistant kernel directionality specified, use resistance not conductance.')
				sys.exit(-1)
			# Get parameter values
			dirtype_A = float(tempdirtype[1])
			dirtype_B = float(tempdirtype[2])
		dirtype = tempdirtype[0]
		# Get scale values for wind and hiking
		if dirtype == 'Hiking' or dirtype == 'Wind':
			if len(minmaxres.split(';')) != 2:
				print('Minimum and maximum scaling values needed for speed to resistance calculation.')
				sys.exit(-1)
			else:
				minres = float(minmaxres.split(';')[0])
				maxres = float(minmaxres.split(';')[1])
		
	# For symmetric models
	if not directionans:
		# Only resistance right now
		if not resans:
			print('Conductance option is not functioning currently with symmetric models. Use resistance.')
			sys.exit(-1)
	
	# ------------------------
	# No direction
	# ------------------------
	# Just resistance surface
	if not directionans:
				
		# Prepare resistance grid file
		resgrid = FileIO.PrepTextFile(gridfile)
		
		# Call the function to convert the raster file to a graph
		#nWeightPairs, xrasloc,yrasloc, header_dict = RasterToNWeights(resgrid)
		nWeightPairs, xrasloc,yrasloc, header_dict = RasterToNWeights2(resgrid)
		
		# Extract information
		nrows = int(header_dict['nrows'])
		ncols = int(header_dict['ncols'])
		xllcorner = float(header_dict['xllcorner'])
		yllcorner = float(header_dict['yllcorner'])
		cellsize = int(header_dict['cellsize'])
		NODATA_value = float(header_dict['NODATA_value'])
	
	# --------------------------
	# Direction functions
	# --------------------------
	else:		
		# Prepare resistance grid files
		resgrid = FileIO.PrepTextFile(gridfile)
		elevgrid = FileIO.PrepTextFile(elevfile) # Or V
		barrgrid = FileIO.PrepTextFile(barrfile) # Or U
		
		if dirtype == 'DEM':
			print('Warning: DEM direction may produce incorrect graph. Recommend using flow accumulation.')
			sys.exit(-1)
		
		# Call the function to convert rasters to a directional graph
		if resans: #Use resistance surface
			if dirtype == 'Wind': # Sum up windspeeds
				nWeightPairs, xrasloc,yrasloc, header_dict = RasterToWindSpeedDirection(resgrid,elevgrid,barrgrid,maxres,minres)
			elif dirtype == 'FlowAcc': # Traditional least cost direction: for riverine and barriers
				nWeightPairs, xrasloc,yrasloc, header_dict = RasterToBarrierDirection(resgrid,elevgrid,barrgrid,dirtype)
			elif dirtype == 'Hiking': # 
				nWeightPairs, xrasloc,yrasloc, header_dict = RasterToHikingWeights(resgrid,elevgrid,maxres,dirtype_A,dirtype_B,minres)
		else: # Use conductance surface
			nWeightPairs, xrasloc,yrasloc, header_dict = RasterToBarrierDirection_Conductance(resgrid,elevgrid,barrgrid,dirtype)
		
		# Extract information
		nrows = int(header_dict['nrows'])
		ncols = int(header_dict['ncols'])
		xllcorner = float(header_dict['xllcorner'])
		yllcorner = float(header_dict['yllcorner'])
		cellsize = int(header_dict['cellsize'])
		NODATA_value = float(header_dict['NODATA_value'])
	
	# Print to log
	start_time = datetime.datetime.now()
	stringout = 'G creation run-time: '+str(datetime.datetime.now() -start_time) + '\n\n'
	FileIO.logMsg(logfHndl,stringout)
	# Timing events: start
	
	# Call the function to read in xy points and get dict. keys for G
	tupMtoRas = MatchToRaster(xyfilename,xrasloc,yrasloc,directionans)
	
	# Retrieve return values here from tuple above
	stringpts = tupMtoRas[0]
	nopoints = tupMtoRas[1]
	
	# Print to log
	start_time = datetime.datetime.now()
	stringout = 'xy points run-time: '+str(datetime.datetime.now()-start_time)+\
	'\n\n'
	FileIO.logMsg(logfHndl,stringout)
	# Timing events: start	

	# Sum of the paths through the resistance grid
	pathadd = np.zeros((nrows,ncols),float)
	if directionans:
		pathadd_rev = np.zeros((nrows,ncols),float)
		pathadd_lower = np.zeros((nrows,ncols),float)
		pathadd_upper = np.zeros((nrows,ncols),float)
	# The cost distance matrix
	cd_matrix = np.zeros((nopoints,nopoints),float)
	if directionans:
		#cd_matrix_rev = np.zeros((nopoints,nopoints),float)
		cd_matrix_rev = np.empty((nopoints,nopoints))
		cd_matrix_rev.fill(np.nan)
		cd_matrix_lower = np.empty((nopoints,nopoints))
		cd_matrix_lower.fill(np.nan)
		cd_matrix_upper = np.empty((nopoints,nopoints))
		cd_matrix_upper.fill(np.nan)
		if dirtype == 'Wind' or dirtype == 'Hiking':
			cd_matrix_lower_accCost = np.empty((nopoints,nopoints))
			cd_matrix_lower_accCost.fill(np.nan)
			cd_matrix_upper_accCost = np.empty((nopoints,nopoints))
			cd_matrix_upper_accCost.fill(np.nan)
	
	# Arrays of the resitance grid x and y values  
	xvalues = np.arange(xllcorner+(cellsize/2.), xllcorner+(cellsize * ncols),cellsize,\
	dtype= 'float')	
	yvalues = np.arange(yllcorner+(cellsize/2.), yllcorner+(cellsize*nrows),cellsize,\
	dtype = 'float')
	
	# yvalues should run from upper to lower to match the r.kwdGetValue('Grid_Filename') and the 
	#	nWeightPairs array, needed for arcgis format
	yvalues = np.flipud(yvalues)
	
	# Index key to track sets of points we've already visited 
	# with the path algorithm
	visited = np.zeros((nopoints,nopoints), dtype = 'b')
	if directionans:
		visited_rev = np.zeros((nopoints,nopoints), dtype = 'b')
		visited_lower = np.zeros((nopoints,nopoints), dtype = 'b')
		visited_upper = np.zeros((nopoints,nopoints), dtype = 'b')
					
	# Store the paths (source, destination, length, paths)
	paths = []
	if directionans:
		paths_rev = []	
		paths_lower = []
		paths_upper = []
	#---------------------------------------------------------------------------
	# End Main loop input variables
	#---------------------------------------------------------------------------      
	#---------------------------------------------------------------------------
	# Start Main Loop Logic
	#---------------------------------------------------------------------------
	
	# Timing events: start
	start_time = datetime.datetime.now()
	# Check for OS environment, if windows run serial_paths:
	if os.name == 'nt' or num_of_pro == 1: 
		stringout = 'Unicor does not currently support multi-processing in the windows os, reverting to serial calculations\n'
		FileIO.logMsg(logfHndl,stringout)
		if not directionans:
			tupSpaths = serial_paths(dirtype,logfHndl,nopoints,\
			stringpts,visited,EDthresholdans,nWeightPairs,\
			xvalues,yvalues,pathadd,cd_matrix,paths,\
			edge_type,edge_dist,transform_func,const_kernal_vol,\
			vol_constant,nbhd_dist,True)
		
		else: # If running direction
			tupSpaths_lower = serial_paths(dirtype,logfHndl,nopoints,stringpts,visited_lower,EDthresholdans,nWeightPairs,xvalues,yvalues,pathadd_lower,cd_matrix_lower,paths_lower,edge_type,edge_dist,transform_func,const_kernal_vol,vol_constant,nbhd_dist,resans,cd_matrix_lower_accCost)
			
			stringpts_rev = stringpts[::-1]
			tupSpaths_upper = serial_paths(dirtype,logfHndl,nopoints,stringpts_rev,visited_upper,EDthresholdans,nWeightPairs,xvalues,yvalues,pathadd_upper,cd_matrix_upper,paths_upper,edge_type,edge_dist,transform_func,const_kernal_vol,vol_constant,nbhd_dist,resans,cd_matrix_upper_accCost)
			
	# If not windows, then run parallel paths
	elif os.name != 'nt' and num_of_pro != 1:
		if directionans:
			tupSpaths_lower = parallel_paths(logfHndl,nopoints,\
			stringpts,visited_lower,EDthresholdans,num_of_pro \
			,nWeightPairs,xvalues,yvalues,pathadd_lower,cd_matrix_lower,paths_lower,edge_type, \
			edge_dist,transform_func,const_kernal_vol, vol_constant,nbhd_dist)

			stringpts_rev = stringpts[::-1]
			tupSpaths_upper = parallel_paths(logfHndl,nopoints,\
			stringpts_rev,visited_upper,EDthresholdans,num_of_pro \
			,nWeightPairs,xvalues,yvalues,pathadd_upper,cd_matrix_upper,paths_upper,edge_type, \
			edge_dist,transform_func,const_kernal_vol, vol_constant,nbhd_dist)
			#print('Parallel processing not implemented with directional cost calculations at this point')
			#sys.exit(-1)
		else:
			tupSpaths = parallel_paths(logfHndl,nopoints,\
			stringpts,visited,EDthresholdans,num_of_pro \
			,nWeightPairs,xvalues,yvalues,pathadd,cd_matrix,paths,edge_type, \
			edge_dist,transform_func,const_kernal_vol, vol_constant,nbhd_dist) 	
	
	# Retrieve return values here from tuple above
	if not directionans:
		pathadd = tupSpaths[0]
		cd_matrix = tupSpaths[1]
		paths = tupSpaths[2]
	else:
		#pathadd_rev = tupSpaths_rev[0]
		#cd_matrix_rev = tupSpaths_rev[1]
		#paths_rev = tupSpaths_rev[2]
		pathadd_lower = tupSpaths_lower[0]
		cd_matrix_lower = tupSpaths_lower[1]
		paths_lower = tupSpaths_lower[2]
		pathadd_upper = tupSpaths_upper[0]
		cd_matrix_upper = tupSpaths_upper[1]
		paths_upper = tupSpaths_upper[2]
		# Lower is ordered 1,2. Upper is reversed 2,1
		# User Lower lower triangle, upper lower triangle to fill upper triangle.
		# Careful, stringpts and stringpts_rev are reversed!
		#if dirtype == 'Wind': # Just for wind also get cost
		#	cd_matrix_upper_accCost = tupSpaths_upper[3]
		#	cd_matrix_lower_accCost = tupSpaths_lower[3]
			
	
	# Print to log
	stringout = '\nTotal shortest path calculation run-time: '+\
	str(datetime.datetime.now() - start_time)+'\n'
	FileIO.logMsg(logfHndl,stringout)
		
	# Get header information, change to -9999 for nodata.
	header_dict, data_list = FileIO.loadFile(resgrid,header_lines=6)
	res_array = np.asarray(data_list,'float')
	pathadd[res_array==-9999] = -9999	
	# --------------------------------------------------------------------------
	# Create output files, and write out important fields
	# --------------------------------------------------------------------------

	# Do some output filename creation
	cd_matrix_ext = '.cdmatrix.csv'    
	pathadd_ext = '.addedpaths.txt'
	paths_ext = '.paths.csv'
	buff_ext = '.kdepaths'
	lev_ext = '.levels'
	Gmets_ext = '.graphmetrics'

	# Strip .rsg and .xy from file names
	input_file_name = gridfile.rstrip('.rsg')
	input_point_name = xyfilename.rstrip('.xy')
	
	# If data is in another directory, clean up more
	if len(input_file_name.split('/')) > 1:
		data_dir = input_file_name.split('/')[0]+'/'
	else:
		data_dir = ''
	input_file_name = input_file_name.split('/')[-1]
	input_point_name = input_point_name.split('/')[-1]	
		
	# Output file names defined here
	if not directionans:
		cd_matrix_file_name = input_file_name + '_' + input_point_name + cd_matrix_ext
		pathadd_file_name = input_file_name  + '_' + input_point_name + pathadd_ext
		paths_file_name = input_file_name + '_' + input_point_name + paths_ext
	else:
		if resans: # Resistance
			if dirtype == 'Wind':
				cd_matrix_file_name = barrfile + '_' + input_file_name + '_' + input_point_name + '_asymmetric_accCost'+cd_matrix_ext
				cd_matrix_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + '_asymmetric_WindCost'+cd_matrix_ext 
				pathadd_file_name_rev_lower = barrfile + '_' + input_file_name  + '_' + input_point_name + '_asymmetric_WindCost_FROM_'+ pathadd_ext
				pathadd_file_name_rev_upper = barrfile + '_' + input_file_name  + '_' + input_point_name + '_asymmetric_WindCost_TO_'+ pathadd_ext				
				paths_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + '_asymmetric_WindCost'+ paths_ext
			else:
				cd_matrix_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + '_asymmetric'+cd_matrix_ext 
				pathadd_file_name_rev = barrfile + '_' + input_file_name  + '_' + input_point_name + '_asymmetric'+ pathadd_ext 
				paths_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + '_asymmetric'+ paths_ext 
		else: # Conductance
			cd_matrix_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + 'asymmetric.condmatrix.csv' 
			pathadd_file_name_rev = barrfile + '_' + input_file_name  + '_' + input_point_name + 'asymmetric'+ pathadd_ext 
			paths_file_name_rev = barrfile + '_' + input_file_name + '_' + input_point_name + 'asymmetric'+ paths_ext
	buff_file_name = input_file_name + '_' + input_point_name + buff_ext
	level_file_name = input_file_name + '_' + input_point_name + lev_ext
	Gmets_file_name = input_file_name + '_' + input_point_name + Gmets_ext

	if DEBUG_MODE:
		#for debugging purposes
		DBug.displayPaths(ncols,nrows,xrasloc,yrasloc,pathadd,resgrid)
	
	# Output cost distance matrix function
	if CDmatrixans: 
		if edge_type.lower() == 'all_paths':
			print('Cannot produce a cost distance matrix for all paths.')
		else: 
			if not directionans:
				output_cdmatrix(logfHndl,cd_matrix,cd_matrix_file_name,data_dir,stringpts)
			else:						
				
				output_cdmatrix_asymmetric(logfHndl,cd_matrix_lower,cd_matrix_upper,cd_matrix_file_name_rev,data_dir,stringpts,resans)
				#if dirtype == 'Wind':
				#	output_cdmatrix_asymmetric(logfHndl,cd_matrix_lower_accCost,cd_matrix_upper_accCost,cd_matrix_file_name,data_dir,stringpts,resans)
				
	
	# Output ascii path addition function
	if outputPathadd:
		if not directionans:
			output_pathadd(logfHndl,pathadd,pathadd_file_name,header_dict,keys,data_dir)
		else:
			if dirtype == 'Wind':
				output_pathadd(logfHndl,pathadd_lower,pathadd_file_name_rev_lower,header_dict,keys,data_dir)
				output_pathadd(logfHndl,pathadd_upper,pathadd_file_name_rev_upper,header_dict,keys,data_dir)
			else:			
				output_pathadd(logfHndl,pathadd_lower + pathadd_upper,pathadd_file_name_rev,header_dict,keys,data_dir)
			
	# This function is to output the paths (node locations) for source
	# 	to destinamtion
	if outputPaths:
		if not directionans:
			output_paths(logfHndl,paths,paths_file_name,data_dir)
		else:
			#output_paths(logfHndl,paths_rev,paths_file_name_rev,data_dir)
			output_paths(logfHndl,paths_lower+paths_upper,paths_file_name_rev,data_dir)		
	
	# This function option outputs the buffered paths (ascii file format)
	if outputBuffer:
		if directionans:
			print('Buffer paths not operational with direction.')
			sys.exit(-1)
		buff_values = createBuffers(ncols,nrows,xrasloc,yrasloc,\
		pathadd,resgrid,BufferGridSize,KernelFunction)
		
		# Output call
		output_buffer(logfHndl,buff_values,buff_file_name,header_dict,keys,data_dir)
	
	# This function option outputs the catagorical buffered paths in 
	# 	ascii file format
	if outputLevels:
		if directionans:
			print('Categorical paths not operational with direction.')
			sys.exit(-1)
		# Case when r.kwdGetValue('Save_KDE_Output') is False
		if not outputBuffer:
			buff_values = createBuffers(ncols,nrows,xrasloc,yrasloc,\
			pathadd,resgrid,BufferGridSize,KernelFunction)
		level_values = createLevels(ncols,nrows,xrasloc,yrasloc,\
		resgrid,buff_values,LevelNumber)
		
		# Output call
		output_levels(logfHndl,level_values,level_file_name,header_dict,keys,data_dir)
	
	# Calculate graph theory metrics
	if outputGraphMetrics:
		
		print("Graph Metrics not operational currently.")
		if False:
		
			# Timing events: start
			start_time = datetime.datetime.now()
			
			# If number of points is greater than 100, don't run
			if nopoints > noptsthreshold:
				
				# Print to log
				stringout = 'Path graph metrics not run do to large graph.\n'
				FileIO.logMsg(logfHndl,stringout)
			
			# If number of points is less than 100, run
			else:
				# Turn -9999 values back to nodata 0 values
				pathadd[pathadd==0.0] = -9999.0
				output_graphmetrics(pathadd,paths,Gmets_file_name,data_dir)
			
				# Print to log
				stringout = 'Path graph metrics completed in: '+str(datetime.datetime.now()-\
				start_time)+'\n'
				FileIO.logMsg(logfHndl,stringout)	
	#-------------------------------------------------------------------------------        
	# END main loop
	#-------------------------------------------------------------------------------

	# Write to log file
	
	FileIO.logMsg(logfHndl,'Total UNICOR program run-time: '+str(datetime.datetime.now()-\
	start_time0),msgVerbose = True)
	logfHndl.close()
	
	return cd_matrix


if __name__ == '__main__':
	
		# --------------------------------------------------------
	# Parse command line, to fetch session input file
	# --------------------------------------------------------
	if len(sys.argv) >= 2:
		ripFilePath = sys.argv[1]
		# If user did not specify .rip file
	else:
		print "User must specify input file name (e.g., at command line type UNICOR.py user_input.rip)."
		sys.exit(-1)
	main(ripFilePath)
