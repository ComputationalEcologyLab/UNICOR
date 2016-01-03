# -------------------------------------------------------------------------------------------------------------
# FileIO.py
# Author(s): Brian Hand
#
# originally created  : May 2011
# GOAL: Basic file I/O that corresponds to the UNICOR and CDPOP standard as of May 2011
# Keywords: File input, output, formatting, ASCII
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
	
# ------------------------------------------------------------------------------------------
#def loadFile(filename, header_lines=0, delimiter=None):
def loadFile(filename, header_lines=0, delimiter=None, cdpop_inputvars=False): ###
	'''
	Used to load file hearders according to current UNICOR standards
	as of 5/19/2011
	'''
	try:
		inputfile = open(filename)
	except (IOError,OSError) as e:
		print("Load file: %s the file (%s) is not available!"%(e,filename))
		sys.exit(-1)
	header_dict = {}
	data_list = []
	index_list = [] ###
	
	if delimiter != None:
		lines = [ln.rstrip().split(delimiter) for ln in inputfile]
	else:
		lines = [ln.rstrip().split() for ln in inputfile]
	
	inputfile.close()
		
	for i,line in enumerate(lines):
		if i < header_lines:
			if len(line) <= 1:
				print("Only one header value in line, skipping line...")
				continue
			#else:	
			elif len(line) == 2: ###
				header_dict[line[0]] = line[1]
			### working with a file where the first line is all header keys
			### and the following lines are their data
			elif cdpop_inputvars: ###
				for j in range(len(line)): ###
					header_dict[line[j]] = [] ###
					index_list.append(line[j]) ###
		else:
			#This is a fix to remove empty entries from from a line if they
			#exist, this is to fix a problem with reading in cdmatrices
			for i in range(line.count('')):
				line.remove('')
			data_list.append(line)
			if cdpop_inputvars: ###
				#tempTuple = ()
				for j in range(len(line)): ###
					if line[j].lower() == 'y':
						line[j] = True
					elif line[j].lower() == 'n':
						line[j] = False
					#---
					# remove the following lines should the tuple representing bar-delimited values break anything -TJJ
					elif line[j].find('|') != -1:
						tempList = line[j].split('|')
						line[j] = tuple(tempList)
					#---
					header_dict[index_list[j]].append(line[j]) ###

	if not cdpop_inputvars:
		return header_dict, data_list
	else:
		n_jobs = len(lines) - header_lines
		return header_dict, index_list, n_jobs
	
# --------------------------------------------------------------------------
def logMsg(outf,msg,msgVerbose = False):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	 msgVerbose -when set True, routes session log traffic to BOTH the
	 screen and to the log file. When False, log traffic just
	sent to log file alone.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	YET-TO-DO: as of 2010-11-17:13:10:01 --refactor logging facility to
	allow working in parallel processing mode.
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print("%s"%(msg))
	
	# End::logMsg()
	
# ------------------------------------------------------------------------------------------

def outputGrid(out_filename,data,headers=None,header_dict=None,keys=None,delimiter=None):
	'''
	Writes to ascii file with header information according to
	current UNICOR standards as of 5/19/2011
	data can either be a python list (or list of lists) or a numpy array (or matrix)
	'''
	########
	if delimiter == None:
		delimiter = ' '
	########
	
	# Create file to write info to
	fout = open(out_filename, 'w')
	if header_dict != None:
		if keys == None:
			keys = header_dict.keys()
		for key in keys:
			#fout.write(key + ' ' + header_dict[key] + '\n')
			fout.write(key + delimiter + header_dict[key] + '\n')
	if headers != None:
		for header in range(len(headers)):
			fout.write(headers[header])
			if header != (len(headers)-1):
				fout.write(delimiter)
		fout.write('\n')
	# Begin loop through rows
	for i in xrange(len(data)):
		# Begin Loop through columns
		for j in xrange(len(data[0])):
			# Write out pathadd list value
			fout.write(str(data[i][j]))
			if j != (len(data[i])-1):
				fout.write(delimiter)
			#if delimiter != None:
			#	fout.write(delimiter)
			#else:
			#	fout.write(' ')
		# Add return character on end
		fout.write('\n')
	# Close file
	fout.close

# ------------------------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files.
	2010-11-10T16:46:01MST jmg --modified to test if textpath
	really a string, as RipMgr kwdGetValue(key) returns native
	Python datatypes and typically not strings...
	'''
	rsltStr = ''
	if type(textpath) == type(str()):

		rsltStr = textpath.rstrip('\n').rstrip('\r')
	else:
		rsltStr = textpath
	return rsltStr 
	# End::PrepTextFile()
	
# --------------------------------------------------------------------------

	def MatchToRaster(xyfilename,xrasloc,yrasloc):
		'''
		MatchToRaster() Match x,y locations to the input grid using this 
		function
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
			s1 = str(float(xrasloc[fixed_pts[1][i]]))
			s2 = str(float(yrasloc[fixed_pts[1][i]]))
			
			# concatenate two additive terms...
			stringpts.append(s1+ '_' + s2)
			# end:: i fixed_pts[] loop
		
		# Return values
		tupMtoG = stringpts, len(stringpts)
		return tupMtoG
		
		# End::MatchToRaster()
# --------------------------------------------------------------------------
#if __name__ == '__main__':


