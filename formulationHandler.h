#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <vector>
#include <list>
#include <array>
#include <Eigen/Dense>
#include <iterator>

#include "gurobi_c++.h"
using namespace std;
using namespace Eigen;

int readSol(string solfilename, int **Xtovec, map<string,int> varNameToNumber, int n, int N, VectorXd *out_vec);
void flushConstraints(GRBModel *mlin, int *M_ptr, int M_base, int *total_cuts_ptr, int N, MatrixXd Abasic, GRBVar *fullx, VectorXd bbasic, tuple<int,string> *const_number);
void boundTightening(GRBModel *m, GRBVar* varlist, int n, map<string,int> varNameToNumber, map<int,string> varNumberToName, int addAll);
GRBModel* linearize(GRBModel *m, map<string,int> varNameToNumber, map<int,string> varNumberToName, bool wRLT, int addAll,
	GRBVar ***out_x, GRBVar ****out_X, bool ***out_isInOriginalModel);
GRBModel* unlinearize(GRBModel *m, GRBVar **x, GRBVar ***X, int n, int M);
void addRLTconstraints(GRBModel *m, GRBVar **x, GRBVar*** X, int n, bool wRLT);
void createMap(GRBVar **x, GRBVar ***X, int ***out_Xtovec, vector< array<int, 2> > *out_vectoX, int n);
void buildAb(GRBModel *m, GRBVar **x, GRBVar ***X, int **Xtovec, int n, vector<RowVectorXd> *out_A, vector<double> *out_b, RowVectorXd *out_c);
void projectDown(GRBModel *m, GRBVar **x, GRBVar ***X, int n, int M, bool **isInOriginalModel, bool keepRLT);
