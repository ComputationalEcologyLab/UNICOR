# -------------------------------------------------------------------------------------------------------------
# RandomFun.py
# Author(s): EL Landguth
#
# originally created  : May 2011
# GOAL: Random functions used as of May 2011
# Keywords: random, weighted draw, probability lists
# ------------------------------------------------------------------------------------------------------
import pdb, random
	
try:
	import numpy as np
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)
# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list. Probilities do not have to add to one.
	List must have 2 entries, an id and a float probability.
	Output returns both the id and the index in the list.
	'''
	wtotal=sum(x[1] for x in lst)
	n=random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	# The case where all of the values in lst are the same
	if len(lst) == count:
		count = count-1
	return item,count
	
	#End::w_choice_general()

def approxWhere(chkAry, searchValue, TOLERANCE = 10E-5):
	''' 
	Uses a where search from numpy on an array to find
	places where a value is located using a 
	tolerance value
	'''
	criterion = (chkAry >= searchValue - TOLERANCE) & \
	(chkAry <= searchValue + TOLERANCE) 
	indices = np.where(criterion)
	return indices

def compareElementsPercent(list1,list2,TOLERANCE = 10E-5):
	'''
	Do an element pair-wise comparision of lists
	using the numpy diff() and where() functions to find
	the differences in each element in the list and return
	a % of the elements that are equal using a tolerance
	'''
	test = abs(np.diff([list1,list2],axis=0))
	if len(np.where(test > TOLERANCE)[0])> 0:
		count_differences = len(np.where(test > TOLERANCE)[0])
	else: count_differences = 0.0
	percent = float(len(list1) - count_differences)/float(len(list1))
	return percent
	
	
if __name__ == '__main__':

	list1 = [2,1,5,6,7]
	list2 = [1,2,3,4,7]
	
	TOLERANCE = 0.25
	percent = compareElementsPercent(list1,list2,TOLERANCE)
	print percent
