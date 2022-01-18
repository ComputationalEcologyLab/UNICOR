# -------------------------------------------------------------------------------------------------------------
# UnicorScript.py
# Author(s): Brian Hand
# originally created  : August 2011
# GOAL: Automate Unicor Runs
# Keywords: Cdmatrix creation, Unicor, automation
# ------------------------------------------------------------------------------------------------------
appRele = ""
appName = "UnicorScript"
appVers = "0.9"
authorNames = "Brian Hand"

#File absolute paths for importing functions
UTILITIES_PATH =  "../utilities/"

# Import Modules with Except/Try statements
# For debugging
import pdb

# Platform and system functions
try:
	import os, sys, random, math
	from operator import attrgetter                    
except ImportError as eMsg:
	print(("ImportError: (%s) OS and SYS and random required!"%(eMsg)))
	import sys
	sys.exit(-1)
	
#Import the package specific folders
CSV_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"CSVParse"))
FILEIO_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"FileIO"))
LOG_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"Log"))
RANDOM_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"RandomFun"))

if CSV_folder not in sys.path:
     sys.path.insert(0, CSV_folder)
if FILEIO_folder not in sys.path:
     sys.path.insert(0, FILEIO_folder)
if LOG_folder not in sys.path:
     sys.path.insert(0,LOG_folder)
if RANDOM_folder not in sys.path:
     sys.path.insert(0,RANDOM_folder)
	 
# Numpy and SciPy Functions
try:
	import numpy as np
	import scipy.stats as spstats
except ImportError as eMsg:
	print(("ImportError (%s) Numpy required."%(eMsg)))
	sys.exit(-1)

# UNICOR functions
try:
	from UNICOR import main as RunUNICOR
except ImportError as eMsg:
	print(("ImportError: (%s) UNICOR.py is required"%(eMsg)))
	sys.exit(-1)

# Utility functions
try: 
	import FileIO
	import CSVParse 
	import Log
	import RandomFun
	from RipMgr import mgrParseStanzaInputs    # Joe's Rip Manager for inputs
	from RipMgr import *
except ImportError as eMsg:
	print(("ImportError: (%s) FileIO.py and CSVParse.py are required."%(eMsg)))
	sys.exit(-1)



def main(*arguements):
	
	# --------------------------------------------------------
	# Parse command line, to fetch session input file
	# --------------------------------------------------------
	if len(sys.argv) >= 2:
		ripFilePath = sys.argv[1]		
	
	# If user did not specify .rip file
	else:
		print("User must specify input file name (e.g., at command line type GARM.py user_input.rip).")
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(ripFilePath):
		print(("Cannot find or open runtime inputs file(%s)"%(ripFilePath)))
		sys.exit(-1)

	dir = 'tree100perlandscape/'

	#bird_types = ['St','Mu','Gl','Cr','Ax',Sh,Tu
	#Skipped Du, Re, Ru
	bird_types = ['Ha']
	file_types = ['_tree_100_1_10.asc','_tree_100_1_100.asc']
	#file_types = ['_tree_100_1_2.asc','_tree_100_1_5.asc','_tree_100_1_10.asc','_tree_100_1_100.asc']
	grid_file_out = 'aus_bird_unicor.rip'
	for bird in bird_types:
		for type in file_types:
			print('now running ' + bird + type)
			header_dict, data = FileIO.loadFile(ripFilePath,header_lines = 16)
			header_dict['XY_Filename'] = str(bird + '_XY.csv')
			header_dict['Grid_Filename'] = str(dir + bird + type)
			header_dict['Session_label'] = bird + '_type'
			FileIO.outputGrid(grid_file_out,['null','null'],header_dict = header_dict) 
			cd_matrix = RunUNICOR(grid_file_out)
if __name__ == '__main__':
	main(*sys.argv[1:])

