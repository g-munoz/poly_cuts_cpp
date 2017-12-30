#include <iostream>
#include <fstream>
#include <map>
#include <cmath> 
#include <vector>
#include <list>
#include <array>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

void generalizedminorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation);

void minorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation);

void outerpsd(VectorXd solX, MatrixXd solX_matrix, int n, int N, int **Xtovec,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation);

void eymcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic,
		RowVectorXd *out_pi, double *out_pirhs, double *out_violation);

void shiftedconeeymcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic,
		RowVectorXd *out_pi, double *out_pirhs, double *out_violation);

void findRays(MatrixXd Abasic, int **Xtovec, int n, int N, vector<VectorXd> *out_dirs, vector<MatrixXd> *out_dirs_matrix);
MatrixXd buildMatrixSol(VectorXd solX, int **Xtovec, int n);
VectorXd buildSolFromMatrix( MatrixXd solX_matrix, int **Xtovec, int n);
int checkifparallel(RowVectorXd pi, double pirhs, vector<RowVectorXd> A, vector<double> b, int M, vector<double> rowNorms);
