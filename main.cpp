#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include <unistd.h>
#include <cstdlib>
#include <Eigen/Dense>
#include "formulationHandler.h"
#include "cutGenerators.h"
#include "gurobi_c++.h"
#include <getopt.h> // getopt, getopt_long

using namespace std;
using namespace Eigen;

#define no_argument 0
#define required_argument 1
#define optional_argument 2

string fullfilename;
string solfilename;

/*misc flags*/
bool doBoundTightening = true;
bool checksol = false;
bool outputLast = false;
bool doflush = true;
int flush_freq = 15;
int print_freq = 5;
int addAll = 2;

/*cut family flags*/
bool wRLT = false;
bool EYM = false;
bool SEYM = true;
bool OA = true;
bool Minor = true;
bool GenMinor = true;
bool Strengthening = true;

/*tolerances*/
double eps_main = 1E-8;
double stall_tol = 1E-5;
double eps_coeff = 1E-9; //minimum absolute value of a coefficient.
double gurobi_tol = 1E-6;
double instability_tol = 1E-3;

/*limits*/
int max_iter_stall = 50;
int max_iter = 100000;
int max_cuts = 20; //max cuts each subroutine will return
double max_run_time = 30;

RowVectorXd obj_vector;

void computeBasis(GRBModel *m, vector<RowVectorXd> A, vector<double> b, GRBVar **x, GRBVar ***X, int n, int ** Xtovec, int N, int M, MatrixXd *out_Abasic, VectorXd *out_bbasic, tuple<int,string> *const_number);
void processArgs(int argc, char* argv[]);

int main(int argc, char *argv[]){

	processArgs(argc, argv);

	if (addAll != 2){
		EYM = false;
		SEYM = false;
		OA = false;
	}

	cout << "=========" << endl;
	cout << "Cut Flags" << endl;
	cout << "EYM\t" << EYM << endl;
	cout << "SEYM\t" << SEYM << endl;
	cout << "OA\t" << OA << endl;
	cout << "Minor\t" << Minor << endl;
	cout << "G-Minor\t" << GenMinor << endl;
	cout << "=========" << endl;

	cout << fullfilename << endl;

	size_t lastindex = fullfilename.find_last_of(".");
	string filename = fullfilename.substr(0, lastindex);

	GRBEnv *env;
	env = new GRBEnv();
	GRBModel *m = new GRBModel(*env, fullfilename);

	int n = m->get(GRB_IntAttr_NumVars);
	int nInt = m->get(GRB_IntAttr_NumIntVars);

	//m->write(fullfilename);
	//exit(1);
	if(nInt > 0){
		printf("INFO:%s has too many variables or contains integer variable\n\n", fullfilename.c_str());
		exit(0);
	}

	GRBVar *varlist = m->getVars();
	GRBVar **x;
	GRBVar ***X;

	bool **isInOriginalModel;

	map<string,int> varNameToNumber;
	map<int,string> varNumberToName;

	for(int i=0; i < n; i++){
		varNumberToName[i] = varlist[i].get(GRB_StringAttr_VarName);
		varNameToNumber[varlist[i].get(GRB_StringAttr_VarName)] = i;
	}

	if(doBoundTightening){
		boundTightening(m, varlist, n, varNameToNumber,varNumberToName, addAll);
		string out_name = filename + "_BT.lp";
		m->write(out_name);
	}

	GRBModel *mlin = linearize(m, varNameToNumber, varNumberToName, wRLT, addAll, &x, &X, &isInOriginalModel);

	//mlin->write("foo.lp");

	int N = mlin->get(GRB_IntAttr_NumVars);
	int M = mlin->get(GRB_IntAttr_NumConstrs);

	int **Xtovec;
	vector< array<int, 2> > vectoX;

	createMap(x, X, &Xtovec, &vectoX, n);
	GRBVar *fullx = new GRBVar[N];

	for(int i=0; i < n; i++)
		fullx[i] = *(x[i]);
	for(int l=0; l < n; l++)
		for(int k=0; k < l+1; k++){
			if(Xtovec[k][l] == -1) continue;
			fullx[Xtovec[k][l]] = *(X[k][l]);
		}

	vector<RowVectorXd> A;
	vector<double> b;
	RowVectorXd c;

	buildAb(mlin, x, X, Xtovec, n, &A, &b, &c);
	obj_vector = c;

	VectorXd truesol;
	if(checksol){
		solfilename = filename + ".res";
		if(readSol(solfilename, Xtovec, varNameToNumber, n, N, &truesol) != 0){
			cout << "Warning, there was an issue reading the solution. Turning off solution check." << endl;
			checksol = false;
		}
		else{
			cout << "Solution read correctly. It has objective value "<< truesol.dot(c) << endl;
		}
	}

	mlin->getEnv().set(GRB_IntParam_OutputFlag,0);
	mlin->getEnv().set(GRB_IntParam_Threads,1);

	mlin->update();
	mlin->optimize();

	if(mlin->get(GRB_IntAttr_Status) != 2){
		cout << "Error, model could not be solved. Gurobi Status" << mlin->get(GRB_IntAttr_Status) << endl;
		return 1;
	}

	int iter_num = 1;
	double old_val = mlin->get(GRB_DoubleAttr_ObjVal);
	double new_val = old_val;
	double RLT_val = old_val;

	int iter_stall = 0;
	double stall_val = old_val;
	/*
	Count for each type of cut, with:
	#Type 0 = EYM
	#Type 1 = ShiftedEYM
	#Type 2 = PSD
	#Type 3 = Minor
	#Type 4 = GeneralizedMinor
	*/

	int *counts = new int[6];
	for (int i = 0; i < 6; i++)
		counts[i] = 0;

	clock_t start_time = clock();

	double gurobi_time = mlin->get(GRB_DoubleAttr_Runtime);
	double cut_time = 0;
	double pre_time = 0;
	double post_time = 0;
	double run_time = 0;
	int total_cuts = 0;

	std::vector<double> rowNorms;
	for(int i=0; i < M; i++){
		rowNorms.push_back(A[i].norm());
	}

	cout << "Linearized model has " << N << " variables" << endl;
	cout << "===================================================" << endl;
	printf("%5s %15s %15s %6s %6s\n", "Iter", "Max Violation", "Objective", "Cuts", "Time");
	cout << "===================================================" << endl;

	GRBModel *mlin_base = new GRBModel(*mlin);
	int M_base = M;
	tuple<int,string> *const_number = new tuple<int,string>[N];

	while(iter_num < max_iter){
		double objval = mlin->get(GRB_DoubleAttr_ObjVal);
		clock_t start_t = clock();

		MatrixXd Abasic;
		VectorXd bbasic;

		computeBasis(mlin, A, b, x, X, n, Xtovec, N, M, &Abasic, &bbasic, const_number);

		VectorXd xbasic(N);
		for (int i = 0; i < N; i++)
			xbasic(i) = fullx[i].get(GRB_DoubleAttr_X);

		/*
		VectorXd xbasic_true = Abasic.colPivHouseholderQr().solve(bbasic);
		VectorXd diff = xbasic - xbasic_true;
		if(diff.norm()/xbasic_true.norm() > eps_main)
			cout << "Warning: Gurobi v Xbasic " << diff.norm()/xbasic_true.norm() << endl;
		*/
		MatrixXd xbasic_matrix = buildMatrixSol(xbasic, Xtovec, n); //This matrix will be [[1 x][x X]]
		/*
		FullPivLU<MatrixXd> lu(xbasic_matrix);
		if(lu.rank() == 1){
			cout << "Current solution of rank 1, terminating algorithm." << endl;
			break;
		}
		*/

		vector<VectorXd> dirs;
		vector<MatrixXd> dirs_matrix;
		findRays(Abasic, Xtovec, n, N, &dirs, &dirs_matrix);

		int pool_size=0;

		pre_time += ( clock() - start_t ) / (double)(CLOCKS_PER_SEC);
		start_t = clock();

		vector<tuple<RowVectorXd, double, double, int>> cut_tuples;

		if(doflush && (iter_num)%flush_freq == 0){
			flushConstraints(mlin, &M, M_base, &total_cuts, N, Abasic, fullx, bbasic, const_number);
			mlin->update();
			A.clear();
			b.clear();
			buildAb(mlin, x, X, Xtovec, n, &A, &b, &c);
			obj_vector = c;
		}
		if(EYM){
			RowVectorXd pi1;
			double pirhs1;
			double violation1;

			eymcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, &pi1, &pirhs1, &violation1);
			cut_tuples.push_back(make_tuple(pi1, pirhs1, violation1, 0));
			pool_size++;

		}
		if(SEYM){
			RowVectorXd pi1;
			double pirhs1;
			double violation1;

			shiftedconeeymcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, Strengthening, &pi1, &pirhs1, &violation1);
			cut_tuples.push_back(make_tuple(pi1, pirhs1, violation1, 1));
			pool_size++;
		}
		if(OA){
			vector<RowVectorXd> pi_all;
			vector<double> pirhs_all;
			vector<double> violation_all;

			outerpsd(xbasic, xbasic_matrix, n, N, Xtovec, &pi_all, &pirhs_all, &violation_all);
			int oa_cuts = pi_all.size();
			for (int i=0; i<oa_cuts; i++)
				cut_tuples.push_back(make_tuple(pi_all[i], pirhs_all[i], violation_all[i], 2));
			pool_size+=oa_cuts;
		}
		if(Minor){
			vector<RowVectorXd> pi_all;
			vector<double> pirhs_all;
			vector<double> violation_all;

			principalminorcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, max_cuts, Strengthening, truesol, checksol,
					&pi_all, &pirhs_all, &violation_all);

			int minor_cuts = pi_all.size();
			for (int i=0; i<minor_cuts; i++)
				cut_tuples.push_back(make_tuple(pi_all[i], pirhs_all[i], violation_all[i], 3));
			pool_size+=minor_cuts;
		}
		if(GenMinor){
			vector<RowVectorXd> pi_all;
			vector<double> pirhs_all;
			vector<double> violation_all;

			generalizedminorcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, max_cuts, Strengthening, truesol, checksol,
					&pi_all, &pirhs_all, &violation_all);
			int minor_cuts = pi_all.size();
			for (int i=0; i<minor_cuts; i++)
				cut_tuples.push_back(make_tuple(pi_all[i], pirhs_all[i], violation_all[i], 4));
			pool_size+=minor_cuts;
		}

		cut_time += ( clock() - start_t ) / (double)(CLOCKS_PER_SEC);
		start_t = clock();

		if(pool_size == 0){
			cout << "\nAlgorithm finished because current pool of cuts does not cut enough" << endl;
			break;
		}

		//Up to here all cuts are in cut_tuples. So now we add them to the model
		int skipped = 0;

		//Here we sort the cuts accoring to violation
		int added_cuts = 0;
		sort(begin(cut_tuples), end(cut_tuples),cut_comparator());

		double max_viol = get<2>(cut_tuples[0]);
		for(int k=0; k < pool_size && added_cuts < max_cuts; k++){

			RowVectorXd curr_pi = get<0>(cut_tuples[k]);
			double curr_pirhs = get<1>(cut_tuples[k]);
			double curr_violation = get<2>(cut_tuples[k]);
			int curr_type = get<3>(cut_tuples[k]);

			if(curr_violation < gurobi_tol) //sorted by violation
				break;

			for(int i=0; i<N; i++)
				if(abs(curr_pi[i]) < eps_coeff) curr_pi[i] = 0; //round small numbers

			/*
			if(checkifparallel(curr_pi, curr_pirhs, A, b, M, rowNorms) == 1){
				skipped ++;
				continue;
			}
			*/

			if(checksol){
				if(truesol.dot(curr_pi) - curr_pirhs > gurobi_tol){
					cout << "Error! Cut cutting given solution by " << (truesol.dot(curr_pi) - curr_pirhs) << ". Type is " << curr_type << endl;
					break;
				}
			}

			added_cuts++;

			A.push_back(curr_pi);
			b.push_back(curr_pirhs);
			rowNorms.push_back(curr_pi.norm());

			GRBLinExpr cut;
			cut.addTerms(curr_pi.data(), fullx, N);
			cut -= curr_pirhs;
			char cut_name[50];
			sprintf(cut_name, "Type%d_%d", curr_type ,counts[curr_type] );
			mlin->addConstr(cut <= 0, cut_name);
			counts[curr_type]++;
			total_cuts++;
			M++;
		}

		run_time =  (clock() - start_time ) / (double)(CLOCKS_PER_SEC);
		post_time +=  ( clock() - start_t ) / (double)(CLOCKS_PER_SEC);

		if(iter_num==1 || iter_num%print_freq == 0)
			printf("%5d %2.12f %2.12f %6d %3.2f\n",iter_num, max_viol , objval, total_cuts, run_time);

		mlin->update();
		mlin->optimize();

		gurobi_time += mlin->get(GRB_DoubleAttr_Runtime);

		if(mlin->get(GRB_IntAttr_Status) != 2){
			cout << "\nError: model could not be solved. Gurobi Status " << mlin->get(GRB_IntAttr_Status) << endl;
			break;
		}

		iter_num += 1;
		new_val = mlin->get(GRB_DoubleAttr_ObjVal);

		if(abs((new_val - old_val)/(old_val + eps_main)) < stall_tol)
			iter_stall ++;
		//if(abs((new_val - stall_val)/(stall_val + eps_main)) < stall_tol)
		//	iter_stall ++;
		else{
			iter_stall = 0;
			stall_val = new_val;
			if((new_val - old_val)/abs(old_val) < -instability_tol && mlin->get(GRB_IntAttr_ModelSense) > 0 ){
				cout << "Warning!! Potential numerical instability" << endl;
				//break;
			}
			else if((new_val - old_val)/abs(old_val) > instability_tol && mlin->get(GRB_IntAttr_ModelSense) < 0){
				cout << "Warning!! Potential numerical instability" << endl;
				//break;
			}
		}

		if(iter_stall > max_iter_stall){
			cout << "\nAlgorithm finished because objective value did not increase in the last " << max_iter_stall << " iterations" << endl;
			break;
		}
		if(run_time > max_run_time){
			cout << "\nMax runtime reached" << endl;
			break;
		}
		if(added_cuts == 0){
			cout << "\nNo more cuts added" << endl;
			break;
		}
		old_val = new_val;
		dirs.clear();
		dirs_matrix.clear();

	}
	printf("INFO:%s,%2.6f,%2.6f,%d,%d,%d,%d,%d,%d,%3.2f,%3.2f,%1.2f,%3.2f,%3.2f,%3.2f\n\n",
					fullfilename.c_str(), RLT_val, old_val, counts[0], counts[1], counts[2], counts[3], counts[4], iter_num, run_time, gurobi_time, gurobi_time/run_time, pre_time, cut_time, post_time );

	if(outputLast){
		cout << "Writing last LPs" << endl;
		MatrixXd Abasic;
		VectorXd bbasic;
		computeBasis(mlin, A, b, x, X, n, Xtovec, N, M, &Abasic, &bbasic, const_number);

		string out_name = filename + "_final.lp";
		mlin->write(out_name);

		GRBModel *mNL_full = unlinearize(mlin, x, X, n, M);
		out_name = filename + "_final_NL.lp";
		mNL_full->write(out_name);

		flushConstraints(mlin, &M, M_base, &total_cuts, N, Abasic, fullx, bbasic, const_number);
		mlin->update();
		out_name = filename + "_final_flushed.lp";
		mlin->write(out_name);

		GRBModel *mNL = unlinearize(mlin, x, X, n, M);
		out_name = filename + "_final_flushed_NL.lp";
		mNL->write(out_name);
	}

	delete[] const_number;
	delete[] fullx;
	vectoX.clear();
	A.clear();
	b.clear();
	//c.clear();

	for(int i = 0; i< n; i++)
		delete[] Xtovec[i];

	delete[] Xtovec;
	delete[] counts;
	delete[] x;
	for(int i = 0; i< n; i++){
		delete[] X[i];
		delete[] isInOriginalModel[i];
	}
	delete[] X;
	delete[] isInOriginalModel;
	delete[] varlist;
	delete m;
	delete mlin;
	delete mlin_base;
	delete env;

  return 0;
}

void computeBasis(GRBModel *m, vector<RowVectorXd> A, vector<double> b, GRBVar **x, GRBVar ***X, int n, int **Xtovec, int N, int M,
			MatrixXd *out_Abasic, VectorXd *out_bbasic, tuple<int, string> *const_number){

	GRBConstr *constrs = m->getConstrs();
	MatrixXd Abasic(N,N);
	VectorXd bbasic(N);

	Abasic.setZero();
	bbasic.setZero();

	int count = 0;

	for(int j=0; j < M; j++){
		if(constrs[j].get(GRB_IntAttr_CBasis) == -1){
			for(int i=0; i<N; i++){
				Abasic(count,i) = A[j](i);
			}
			bbasic(count) = b[j];
			const_number[count] = make_tuple(j,constrs[j].get(GRB_StringAttr_ConstrName));
			count++;
		}
	}

	if(count < N){
		for(int i=0; i < n; i++){
			if(x[i]->get(GRB_IntAttr_VBasis) == -1){
				Abasic(count,i) = -1;
				bbasic(count) = -x[i]->get(GRB_DoubleAttr_LB);
				const_number[count] = make_tuple(-1, "bound"); //to flag bound
				count++;
			}
			else if( x[i]->get(GRB_IntAttr_VBasis) == -2){
				Abasic(count,i) = 1;
				bbasic(count) = x[i]->get(GRB_DoubleAttr_UB);
				const_number[count] = make_tuple(-1, "bound");
				count++;
			}

			if(count == N) break;
		}

		for(int l=0; l<n; l++){
			for(int k=0; k<l+1; k++){
				if(X[k][l] == nullptr) continue;

				if(X[k][l]->get(GRB_IntAttr_VBasis) == -1){
					Abasic(count, Xtovec[k][l]) = -1;
					bbasic(count) = -X[k][l]->get(GRB_DoubleAttr_LB);
					const_number[count] = make_tuple(-1, "bound");
					count++;
				}
				else if(X[k][l]->get(GRB_IntAttr_VBasis) == -2){
					Abasic(count, Xtovec[k][l]) = 1;
					bbasic(count) = X[k][l]->get(GRB_DoubleAttr_UB);
					const_number[count] = make_tuple(-1, "bound");
					count++;
				}
				if(count == N) break;
			}
			if(count == N) break;
		}
	}

	*out_Abasic = Abasic;
	*out_bbasic = bbasic;

	delete[] constrs;
}

void processArgs(int argc, char* argv[]) {
	const char* const short_opts = "e:s:o:m:f:g:h:b:w:l:p:t:v:a:";
	const struct option long_opts[] =
  		{
		{"eym", 			required_argument,	0,		'e'},
		{"seym", 			required_argument,	0,		's'},
		{"oa", 				required_argument,	0,		'o'},
		{"minor",      		required_argument,  0,      'm'},
    	{"filename",      	required_argument,	0,      'f'},
    	{"genminor",      	required_argument,	0,      'g'},
		{"flush",      		required_argument,	0,      'h'},
		{"boundtightening",	required_argument,	0,      'b'},
    	{"wrlt",     		required_argument,	0,      'w'},
    	{"timelimit",     	required_argument,	0,      'l'},
    	{"outputlast",  	required_argument,	0,		'p'},
    	{"strengthening",  	required_argument,	0,		't'},
    	{"verify",     		required_argument,	0,		'v'},
		{"addall",     		required_argument,	0,		'a'},
    	{nullptr,     		no_argument,        nullptr,  0}
  };

	std::string helpString = "usage: " + std::string(argv[0]) + " -f <filename> [-b 0 | -b 1] [-l timelimit]\n";
	if ( argc < 2 ){
		cout << helpString.c_str();
		exit(1);
	}
	int inp;
	int option_index;
	while ((inp = getopt_long(argc, argv, short_opts, long_opts, &option_index)) != -1) {
			switch(inp){
				case 's':
					if(optarg && atoi(optarg) != 0)	SEYM = true;
					else SEYM = false;
					break;
				case 'e':
					if(optarg && atoi(optarg) != 0)	EYM = true;
					else EYM = false;
					break;
				case 'o':
					if(optarg && atoi(optarg) != 0)	OA = true;
					else OA = false;
					break;
				case 'm':
					if(optarg && atoi(optarg) != 0)	Minor = true;
					else Minor = false;
					break;
				case 'g':
					if(optarg && atoi(optarg) != 0)	GenMinor = true;
					else GenMinor = false;
					break;
				case 'b':
					if(optarg && atoi(optarg) != 0)	doBoundTightening = true;
					else doBoundTightening = false;
					break;
				case 'w':
					if(optarg && atoi(optarg) != 0)	wRLT = true;
					else wRLT = false;
					break;
				case 'h':
					if(optarg && atoi(optarg) > 0){
						doflush = true;
						flush_freq = atoi(optarg);
					}
					else doflush = false;
					break;
				case 'f':
					fullfilename = optarg;
					break;
				case 'p':
					if(optarg && atoi(optarg) != 0)	outputLast = true;
					else outputLast = false;
					break;
				case 't':
					if(optarg && atoi(optarg) != 0)	Strengthening = true;
					else Strengthening = false;
					break;
				case 'v':
					if(optarg && atoi(optarg) != 0)	checksol = true;
					break;
				case 'a':
					if(optarg) addAll = atoi(optarg);
					break;
				case 'l':
					if(optarg && atof(optarg) > eps_main) max_run_time = atof(optarg);
					break;

			}//switch
	} // process arguments
	for (int i = 0; i < argc; i++) {
		cout << argv[i] << " ";
	}
	cout << endl;
} /* processArgs */
