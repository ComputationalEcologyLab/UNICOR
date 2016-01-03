#!/usr/bin/env python
# file: csvParse.py

import sys
import os.path
from numpy import zeros

# check the number of command line arguments
#if len(sys.argv) != 3:
#	print "Incorrect usage; you must supply a file path and file name"
#	print "Usage: csvParse.py <filePath> <csvFile>"
#	sys.exit(-1)
	
#filePath = sys.argv[1]
#csvFile = sys.argv[2]

#================================================================
# csvParse(filePath, csvFile):
#----------------------------------------------------------------
# given a file path and a csv file, returns a matrix of values
#----------------------------------------------------------------
# input:
#       filePath: string literal of directory in which the csv resides
#       csvFile:  string literal of the file to be parsed
# output:
#       a matrix of values parsed from the csv
#================================================================
def csvParse(filePath, csvFile):
  # ensure the file path exists
  if not os.path.isdir(filePath):
    print "Path does not exist!"
    sys.exit(-1)

  # ensure the csv file exists
  if not os.path.isfile(filePath + csvFile):
    print "File does not exist!"
    sys.exit(-1)

  # read the file...
  infile = filePath + csvFile
  file = open(infile,"r")
  lines = file.readlines()
  file.close()

  # determine the number of rows and columns of the matrix
  nrows = len(lines)
  temp_splitline = lines[0].split(',')
  temp_splitline.remove('\n')
  ncols = len(temp_splitline)

  # create an empty matrix[nrows][ncols]
  a = zeros((nrows,ncols))

  for i in range(0,len(lines)):
    line = lines[i]
    splitline = line.split(',')
    splitline.remove('\n')
    for j in range(0,len(splitline)):
      a[i][j] = splitline[j]

  return a

#================================================================
# generatezt(inMatrix,outFile):
#----------------------------------------------------------------
# given a matrix, generate a file that is readable by zt
#----------------------------------------------------------------
# input:
#       inMatrix: the matrix to be converted to the zt format
#       outFile:  sring literal to name the output file
# output:
#       a .zt file (plain-text) to be read by zt, named <outfile>.zt
#================================================================
def generatezt(inMatrix,outFile):
  # determine dimensions of inMatrix
  matDimension = inMatrix.shape
  matDimension = matDimension[0]  # input is a square matrix, so [0]==[1]
  # open file for writing  
  file = open(outFile, 'w')
  if file:
    # write the inMatrix dimensions to the file
    file.write(str(matDimension))
    file.write("\n")
    # grab all values of lower tri and write to file
    index = 0
    for i in range(1,matDimension):
      for j in range(0,i):
        #outMatrix[i][j] = inMatrix[i][j]
        file.write(str(inMatrix[i][j]))
        if not(j==(i-1)):  
          file.write(" ")
        index = index+1
      file.write("\n")
    file.close()
  else:
    print "Error opening output file!"

#================================================================
# csv2zt(filePath,csvFile)
#----------------------------------------------------------------
# given a filePath and a csvFile; generate a .zt-readable file
#----------------------------------------------------------------
#================================================================
def csv2zt(filePath,csvFile):
  matrix = csvParse(filePath,csvFile)
  output = filePath + csvFile.rstrip(".csv") + ".zt"
  generatezt(matrix,output)
  