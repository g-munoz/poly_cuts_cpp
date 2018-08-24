#!/usr/bin/python

from gurobipy import *
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix
import math
import sys
import os
import time

def main():
	
	print ("Reading solution file: %s" % sys.argv[2] )

	inputVals = {}
	with open(sys.argv[2]) as f:
		for line in f:
			lineVal = line.split()
			if lineVal[0][0] == 'x':
				inputVals[lineVal[0]] = float(lineVal[2])

	print ("Loading model file: %s" % sys.argv[1] )
	m = read(sys.argv[1])
	varlist = m.getVars()
	varNameToNumber = {}
	varNumberToName = {}
	for i in xrange(m.numVars):
		name = varlist[i].getAttr("VarName")
		if name[0] == 'X':
			parts = name.split(",")
			name1 = parts[0][2:]
			name2 = parts[1][:len(parts[1])-1]
			m.addConstr(varlist[i] <= (1+eps)*inputVals[name1]*inputVals[name2])
			m.addConstr(varlist[i] >= (1-eps)*inputVals[name1]*inputVals[name2])

			print 'Added', name, '=', inputVals[name1]*inputVals[name2]
		else:
			m.addConstr(varlist[i] <= (1+eps)*inputVals[name])
			m.addConstr(varlist[i] >= (1-eps)*inputVals[name])
			print 'Added', name, '=', inputVals[name]

	m.update()
	filename, file_extension = os.path.splitext(sys.argv[1])
	
	m.optimize()

	if not m.getAttr("Status") == 2:
		print 'Error in', filename
		#m.write("%s_error.lp" % filename)
		print 'Error, model could not be solved. Gurobi Status', m.getAttr("Status")
		
		sys.exit(2)
	
	
if __name__ == "__main__":
	global eps
	eps = 0
	
	main()

