# Copyright 2017, Gurobi Optimization, Inc.

GUROBI_DIR = /opt/gurobi/gurobi752/linux64
INC      = $(GUROBI_DIR)/include/
EIG_LIB  = .
CPP      = g++
CARGS    = -m64 -g -std=c++11 -O3
CPPLIB   = -L$(GUROBI_DIR)/lib/ -lgurobi_c++ -lgurobi75

main: main.o formulationHandler.o cutGenerators.o
	$(CPP) $(CARGS) -o $@ $^ -I$(INC) -I$(EIG_LIB) $(CPPLIB) -lm

cutGenerators.o: cutGenerators.cpp
	$(CPP) $(CARGS) -c $< -I$(INC) -I$(EIG_LIB) $(CPPLIB) -lm 

formulationHandler.o: formulationHandler.cpp
	$(CPP) $(CARGS) -c $< -I$(INC) -I$(EIG_LIB) $(CPPLIB) -lm 

main.o: main.cpp
	$(CPP) $(CARGS) -c $< -I$(INC) -I$(EIG_LIB) $(CPPLIB) -lm


clean:
	rm -rf *.o
