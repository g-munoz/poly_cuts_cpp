#!/usr/bin/python

from pyscipopt import *

import sys
import os

def main():

	m = Model()
	
	x = [None for i in xrange(5)]

	for i in xrange(5):
		x[i] = m.addVar(lb=0, ub=1, name="x(%d)"%(i+1))
	
	objNL = m.addVar(lb=-m.infinity(), name="obj")

	m.setObjective( 42*x[0] + 44*x[1] + 45*x[2] + 47*x[3] + 47.5*x[4] + objNL, "minimize")
	m.addCons(-50*x[0]*x[0] -50*x[1]*x[1] -50*x[2]*x[2] -50*x[3]*x[3] -50*x[4]*x[4] - objNL <= 0)
	m.addCons(20*x[0] + 12*x[1] + 11*x[2] + 7*x[3] + 4*x[4] <= 40)


	## Constraints below are redundant. If all check are marked as "True", the optimal value is the same as the
	## value for the problem defined above. This optimal value is -17
	## However, if any of them is marked as "False" the optimal value jumps to 0.

	checkFlag1 = False
	checkFlag = False

	m.addCons( - 0.1269800912545133*x[0] + 0.1383331882867284*x[1] \
		   + 0.0881458122112916*x[2] + 0.0505088146549664*x[3] \
		   + 0.0721460556206724*x[4] + 0.3720198720697269*x[0]*x[0] \
		   + 2.57326e-05*x[0]*x[1] + 0.0454443963366067*x[1]*x[1] \
		   + 0.0650104900932718*x[2]*x[2] + 0.0413855468527463*x[3]*x[3] \
		   <= 0.5433860511151614, check=checkFlag1)

	m.addCons( - 0.1269795086893034*x[0] + 0.1383325536352482*x[1] \
		   + 0.088142374163711*x[2] + 0.0505085829282609*x[3] \
		   + 0.0721457246256193*x[4] + 0.3720181652995336*x[0]*x[0] \
		   + 2.57325e-05*x[0]*x[1] + 0.0454441878446767*x[1]*x[1] \
		   + 0.0650178133300319*x[2]*x[2] + 0.0413853569821961*x[3]*x[3] \
		   <= 0.5433835581432043, check=checkFlag)

 	m.addCons( - 0.1269313752495149*x[0] + 0.1383313773890496*x[1] \
		   + 0.0881481341835282*x[2] + 0.0505102026019784*x[3] \
		   + 0.0721480381439806*x[4] + 0.3720300949164391*x[0]*x[0] \
		   - 2.6472e-05*x[0]*x[1] + 0.0454750931100892*x[1]*x[1] \
		   + 0.0650125283039153*x[2]*x[2] + 0.0413866840987198*x[3]*x[3] \
		   <= 0.5434248187809282, check=checkFlag)

 	m.addCons( - 0.1272396650689794*x[0] + 0.1382887698068788*x[1] \
		   + 0.0881047991680785*x[2] + 0.0504865414804913*x[3] \
		   + 0.0721100495394896*x[4] + 0.3720127988797262*x[0]*x[0] \
		   + 3.83341e-06*x[0]*x[1] + 0.045413521874779*x[1]*x[1] \
		   + 0.0649769961221167*x[2]*x[2] + 0.0413630246499493*x[3]*x[3] \
		   <= 0.542966113820981, check=checkFlag)

	m.optimize()
	
if __name__ == "__main__":
	
	main()



