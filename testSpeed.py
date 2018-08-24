#!/usr/bin/python

from pyscipopt import *
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix
import math
import sys
import os
import time
import subprocess

def main():

	time_limit = 20.0
	var_limit = 70

	filename, file_extension = os.path.splitext(sys.argv[1])
	times = [0.0, 0.0, 0.0, 0.0, 0.0]

	curr_file = filename+file_extension
	print ("Loading model file: %s" % curr_file )
	m = Model()
	m.readProblem(curr_file)
	m.setRealParam("limits/time", time_limit)

	if len(m.getVars()) >= var_limit:
		print "INFO: Instance",curr_file, "has too many variables:",len(m.getVars()) 
		exit()

	m.optimize()
	if not m.getStatus() == "timelimit":
		print "INFO: SCIP solved",curr_file, "quickly"
		exit()

	
		
	SCIPgap = m.getGap()

	check_file = filename + "_projected_NL" + file_extension
	if not os.path.exists(check_file):
		print "INFO: File", check_file,"does not exist, running code!"
		process = subprocess.Popen("./main -f %s -g1 -h15 -l1"%curr_file, shell=True)
		process.wait()

	curr_file = filename + "_projected_NL" + file_extension
	m.readProblem(curr_file)
	m.setRealParam("limits/time", time_limit)
	m.optimize()
	
	print "INFO:", filename, "SCIP Gap", SCIPgap, "My Gap", m.getGap()
	
if __name__ == "__main__":
	
	main()



