import sys
from datetime import *

def getTime():
	return datetime.now()

class Timer:
	"""Timer: A simple stopwatch for timing processes."""
	def __init__(self):
		self.startTime = datetime.now()
		self.stopTime = self.startTime
		self.delta = 0.0

	def start(self):
		self.startTime = datetime.now()

	def stop(self):
		self.stopTime = datetime.now()
		self.__computeDelta()

	def __computeDelta(self):
		self.delta = self.stopTime - self.startTime

	def getTime(self):
		return self.delta

	def printTime(self):
		print "Time between start and stop: " + str(self.delta)


class Log:
	"""Log: A simple logging utility that incorporates timers."""
	def __init__(self,outFile='default.log',verbose=True,timeStamp=True):
		self.outFile = outFile			#str:the name of output file to be created
		self.verbose = verbose			#bool:whether log messages are sent to the console
		self.timeStamp = timeStamp		#bool:whether log messages are time stamped
		self.__file = ''				#file handle:private, simply the handle to the log file
		self.masterTimer = Timer()		#Timer:tracks entire time a file is open
		self.processTimer = Timer()		#Timer:used for timing blocks of code
		self.__start = True				#bool:private, whether the process timer needs to be started or stopped

	#--------------------------------
	#private, simply returns the time so log entries can be stamped
	#input: none
	#output: the current datetime
	#--------------------------------
	def __timeStamp(self):
		return datetime.now()

	#--------------------------------
	#writes a message to the log
	#input:
	#		msg -- the string to write to log
	#		verboseOverride -- overrides the log's verbose attribute, can only be used to print to console, not to hide a message from console (assuming the attribute is set to False)
	#		timeStampOverride -- same as above, but with regard to time stamps
	#output:
	#		none. message written to log (and possibly console, potentially timestamped)
	#--------------------------------
	def write(self,msg,verboseOverride=False,timeStampOverride=False):
		ts = self.__timeStamp()
		#when verbose, timestamp is ALWAYS written to console, even if not written to log
		if self.verbose or verboseOverride:
			print("[%s] %s"%(ts,msg))
		if self.timeStamp or timeStampOverride:
			print >> self.__file, "[%s]" %(ts),
		print >> self.__file, "%s" % msg

	#--------------------------------
	#writes a message to the log as well as starting/stopping the process timer (for timing/profiling blocks of code)
	#used in pairs... the first call starts timer, second call stops timer... and so on
	#input:
	#		msg -- the string to write to log
	#		verboseOverride -- overrides the log's verbose attribute, can only be used to print to console, not to hide a message from console (assuming the attribute is set to False)
	#		timeStampOverride -- same as above, but with regard to time stamps
	#output:
	#		none. message written to log (and possibly console, potentially timestamped)
	#		elapsed time is also logged, assuming the process timer is stopped.
	#--------------------------------
	def profile(self,msg,verboseOverride=False,timeStampOverride=True):
		if self.__start:
			self.processTimer.start()
			self.__start = False
			#self.write("------------------------------------------------",verboseOverride)
			#write message to log
			self.write(msg,verboseOverride,timeStampOverride)
		else:
			self.processTimer.stop()
			self.__start = True
			#get time delta
			processTime = self.processTimer.getTime()
			#write message to log
			self.write(msg,verboseOverride,timeStampOverride)
			#write time delta to log
			self.write("Elapsed process time: %s" %(processTime),verboseOverride,timeStampOverride)
			#self.write("------------------------------------------------",verboseOverride)

	#--------------------------------
	#opens the log file for writing and sets file handle, or fails and exits
	#input: none
	#output: none
	#--------------------------------
	def open(self):
		#attempt to open outFile for writing
		try:
			self.__file = open(self.outFile,'w')
		except IOerr:
			print("Error opening %s for writing."%(self.outFile))
			sys.exit(-1)
		#start master timer
		self.masterTimer.start()
		#broadcast message to console, timestamp file, write message to file
		self.write("%s opened for writing."%(self.outFile),True,True)
		self.write("================================================",True)

	#--------------------------------
	#closes the file from further writing
	#input: none
	#output: none
	#--------------------------------
	def close(self):
		#stop master timer
		self.masterTimer.stop()
		self.write("================================================",True)
		#show total time file was open (delta master time)
		masterDelta = self.masterTimer.getTime()
		self.write("%s was open for %s"%(self.outFile,masterDelta),True,True)
		#broadcast message to console, timestamp file, write message to file
		self.write("%s closed for writing."%(self.outFile),True,True)
		#close the file from further writing
		self.__file.close()