#!/usr/bin/python

from pyscipopt import *
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix
import math
import sys
import os
import time

def main():

	time_limit = 60.0

	filename, file_extension = os.path.splitext(sys.argv[1])
	times = [0.0, 0.0, 0.0, 0.0, 0.0]

	check_file = filename + "_projected_NL" + file_extension
	if not os.path.exists(check_file):
		print "INFO: File", check_file,"does not exist"	
		exit(2)

	curr_file = filename+file_extension
	print ("Loading model file: %s" % curr_file )
	m = Model()
	m.readProblem(curr_file)
	
	m.setRealParam("limits/time", time_limit)
	start = time.time()
	m.optimize()
	end = time.time()
	times[0] = (end - start)

	curr_file = filename + "_projected_NL" + file_extension
	m.readProblem(curr_file)
	m.setRealParam("limits/time", time_limit)
	start = time.time()
	m.optimize()
	end = time.time()
	times[1] = (end - start)

	if(1):
		m.readProblem(curr_file)
		cons = m.getConss()
		for i in xrange(len(cons)):
			#print cons[i], cons[i].isChecked()
			#print cons[i], cons[i].isEnforced()
			#print cons[i], cons[i].isInitial()
			#print cons[i], cons[i].isRemovable()
			#print SCIPconsIsChecked(cons[i])
			if cons[i].name[:4] == "Type":
				m.delCons(cons[i])
				#m.chgEnforced(cons[i],False)
				#m.chgInitial(cons[i],False)
				#m.chgRemovable(cons[i],True)
				#do1 = False
		#cons = m.getConss()	
		#for i in xrange(len(cons)):
		#	print cons[i], cons[i].isChecked()
		#	print cons[i], cons[i].isEnforced()
		#	print cons[i], cons[i].isInitial()
		#	print cons[i], cons[i].isRemovable()
		#raw_input()
		m.setRealParam("limits/time", time_limit)
		start = time.time()
		m.optimize()
		end = time.time()
		times[2] = (end - start)
		#m.printBestSol()

		if(0):
			curr_file = filename + "_final_flushed_NL" + file_extension
			m.readProblem(curr_file)
			m.setRealParam("limits/time", time_limit)
			start = time.time()
			m.optimize()
			end = time.time()
			times[3] = (end - start)

			m.readProblem(curr_file)
			cons = m.getConss()
			for i in xrange(len(cons)):
				if cons[i].name[:4] == "Type":
					m.chgCheck(cons[i],False)
					m.chgEnforced(cons[i],False)
					m.chgInitial(cons[i],False)
					m.chgRemovable(cons[i],True)

			m.setRealParam("limits/time", time_limit)
			start = time.time()
			m.optimize()
			end = time.time()
			times[4] = (end - start)

	print "INFO", filename, "Orig", times[0], "Cuts", times[1], "Obj", times[2]
	#print "INFO", filename, "Orig", times[0], "Cuts", times[1], "Cuts-NoCheck", times[2], "Flushed", times[3], "Flushed-NoCheck", times[4]
	#print "INFO", filename, "Orig", times[0], "Cuts", times[1]
	
if __name__ == "__main__":
	
	main()



