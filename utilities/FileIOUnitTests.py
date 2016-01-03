# -------------------------------------------------------------------------------------------------------------
# GARMUnitTests.py
# Author(s): Brian Hand
#    
# originally created  : May 2011
# GOAL: Basic unit testing of I/O functions as of May 2011
# Keywords: Unit tests, FileIO
# ------------------------------------------------------------------------------------------------------

from FileIO import *
import unittest

test_file_load = 'test_file.txt'
test_file_save = 'test_file_out.txt'
header_lines = 6

class TestFileIO(unittest.TestCase):
	def setUp(self):
		self.header_dict_out, self.data_out = loadFile(test_file_load, header_lines)
		self.header_dict_in = None
		self.data_in = None
		
	def tearDown(self):
		del(self.header_dict_out)
		del(self.header_dict_in)
		del(self.data_out)
		del(self.data_in)
		
	def testReadWrite(self):
		outputGrid(test_file_save,self.header_dict_out,self.data_out)
		self.header_dict_in, self.data_in = loadFile(test_file_save, header_lines)
		
		self.assertEqual(self.header_dict_in,self.header_dict_out)
		self.assertEqual(self.data_in,self.data_out)
	
if __name__ == '__main__':
    unittest.main()
	