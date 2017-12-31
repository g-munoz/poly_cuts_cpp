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

using namespace std;
using namespace Eigen;


bool doBoundTightening = true;
bool EYM = false;
bool SEYM = true;
bool OA = true;
bool Minor = true;
bool GenMinor = false;
double eps_main = 1E-8;
double stall_tol = 1E-5;

bool doflush = true;
int flush_freq = 10;

void computeBasis(GRBModel *m, vector<RowVectorXd> A, vector<double> b, GRBVar *x, GRBVar **X, int n, int ** Xtovec, int N, int M, MatrixXd *out_Abasic, VectorXd *out_bbasic, int *const_number);

int main(int argc, char *argv[]){
	if ( argc < 2 ){
        	cout<<"usage: "<< argv[0] <<" <filename>\n";
        	return 1;
	}
	int inp;
	string fullfilename;
	while( ( inp = getopt (argc, argv, "e:s:o:m:f:g:h:b:") ) != -1 ){
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
		}
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

	GRBVar* varlist = m->getVars();
	GRBVar *x;
	GRBVar **X;

	map<string,int> varNameToNumber;
	map<int,string> varNumberToName;

	for(int i=0; i < n; i++){
		varNumberToName[i] = varlist[i].get(GRB_StringAttr_VarName);
		varNameToNumber[varlist[i].get(GRB_StringAttr_VarName)] = i;
	}
		
	if(doBoundTightening)
		boundTightening(m, varlist, n, varNameToNumber,varNumberToName);
	
	GRBModel *mlin = linearize(m, varNameToNumber, varNumberToName, false, &x, &X);
	int N = mlin->get(GRB_IntAttr_NumVars);
	int M = mlin->get(GRB_IntAttr_NumConstrs);

	int **Xtovec;
	vector< array<int, 2> > vectoX;
	
	createMap(x, X, &Xtovec, &vectoX, n);
	
	GRBVar *fullx = new GRBVar[N];
	
	for(int i=0; i < n; i++)
		fullx[i] = x[i];
	for(int l=0; l < n; l++)
		for(int k=0; k < l+1; k++)
			fullx[Xtovec[k][l]] = X[k][l];
	
	vector<RowVectorXd> A;
	vector<double> b;
	vector<double> c;

	buildAb(mlin, x, X, Xtovec, n, &A, &b, &c);

	mlin->getEnv().set(GRB_IntParam_OutputFlag,0);

	mlin->update();
	mlin->optimize();
		
	double gurobi_tol = 1E-6;

	if(mlin->get(GRB_IntAttr_Status) != 2){
		cout << "Error, model could not be solved. Gurobi Status" << mlin->get(GRB_IntAttr_Status) << endl;
		return 1;
	}

	int iter_num = 1;
	double old_val = mlin->get(GRB_DoubleAttr_ObjVal);
	double new_val = old_val;

	double RLT_val = old_val;
	
	int iter_stall = 0;
	int max_iter_stall = 10;
	int max_iter = 1000000;

	int *counts = new int[6];
	for (int i = 0; i < 6; i++)
		counts[i] = 0;

	/*
	Count for each type of cut, with:
	#Type 0 = EYM
	#Type 1 = ShiftedEYM
	#Type 2 = PSD
	#Type 3 = Minor
	#Type 4 = GeneralizedMinor
	*/
	
	clock_t start_time = clock();
	double max_run_time = 3600;
	
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
	int *const_number = new int[N];

	//mlin_base->write("base.lp");

	while(iter_num < max_iter){
		double objval = mlin->get(GRB_DoubleAttr_ObjVal);
		clock_t start_t = clock();
		
		MatrixXd Abasic;
		VectorXd bbasic;

		computeBasis(mlin, A, b, x, X, n, Xtovec, N, M, &Abasic, &bbasic, const_number);

		VectorXd xgurobi(N);
		for (int i = 0; i < N; i++){
			xgurobi(i) = fullx[i].get(GRB_DoubleAttr_X);
		}

		VectorXd xbasic = Abasic.colPivHouseholderQr().solve(bbasic);

		VectorXd diff = xbasic - xgurobi;
		if(diff.norm()/xgurobi.norm() > eps_main)
			cout << "Warning: Gurobi v Xbasic " << diff.norm()/xgurobi.norm() << endl;


		//MatrixXd xgurobi_matrix = buildMatrixSol(xgurobi, Xtovec, n); //This matrix will be [[1 x][x X]]
		MatrixXd xbasic_matrix = buildMatrixSol(xbasic, Xtovec, n); //This matrix will be [[1 x][x X]]

		FullPivLU<MatrixXd> lu(xbasic_matrix);
		if(lu.rank() == 1){
			cout << "Current solution of rank 1, terminating algorithm." << endl;
			break;
		}
		
		vector<VectorXd> dirs;
		vector<MatrixXd> dirs_matrix;
		findRays(Abasic, Xtovec, n, N, &dirs, &dirs_matrix);

		int pool_size=0;

		pre_time += ( clock() - start_t ) / (double)(CLOCKS_PER_SEC);
		start_t = clock();

		vector<tuple<RowVectorXd, double, double, int>> cut_tuples;
		//int max_cuts = std::numeric_limits<int>::max();
		int max_cuts = 100; //max cuts each subroutine will return

		if(doflush && (iter_num)%flush_freq == 0){
			GRBConstr *constrs = mlin->getConstrs();
			for(int l=M_base; l < M; l++)
				mlin->remove(constrs[l]);

			mlin->update();
			
			M = M_base;
			total_cuts = 0;

			for(int l=0; l < N; l++)
				if(const_number[l] >= M_base){ //if the constraint is not in the base linearization
					GRBLinExpr basis_row;
					for(int idx = 0; idx < N ; idx++)
						basis_row += Abasic(l,idx)*fullx[idx];
					basis_row -= bbasic[l];
					mlin->addConstr(basis_row <= 0);
					total_cuts++;
					M++;
				}

			mlin->update();
			A.clear();
			b.clear();
			c.clear();
			buildAb(mlin, x, X, Xtovec, n, &A, &b, &c);
			delete[] constrs;
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
			
			shiftedconeeymcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, &pi1, &pirhs1, &violation1);
		
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

			minorcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, max_cuts,
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

			generalizedminorcut(xbasic, xbasic_matrix, dirs, dirs_matrix, n, N, Xtovec, Abasic, bbasic, max_cuts,
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
		
		//Here I'll sort the cuts accoring to violation
		max_cuts = 20; //max cuts added per iteration
		int added_cuts = 0;
		
		sort(begin(cut_tuples), end(cut_tuples),
			  [](tuple<RowVectorXd, double, double, int> const &t1, tuple<RowVectorXd, double, double, int> const &t2) 
			     {return get<2>(t1) > get<2>(t2);} );

		double max_viol = get<2>(cut_tuples[0]);
		
		for(int k=0; k < pool_size && added_cuts < max_cuts; k++){
			
			RowVectorXd curr_pi = get<0>(cut_tuples[k]);
			double curr_pirhs = get<1>(cut_tuples[k]);
			double curr_violation = get<2>(cut_tuples[k]);
			int curr_type = get<3>(cut_tuples[k]);
			
			if(curr_violation < gurobi_tol) //sorted by violation
				break;
			
			/*			
			if(checkifparallel(curr_pi, curr_pirhs, A, b, M, rowNorms) == 1){
				skipped ++;
				continue;
			}
			*/

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
		
		if(iter_num==1 || iter_num%5 == 0)
			printf("%5d %2.12f %2.12f %6d %3.2f \t %6d parallel cuts skipped\n",iter_num, max_viol , objval, total_cuts, run_time, skipped);
		
		
		mlin->update();
		mlin->optimize();
		
		gurobi_time += mlin->get(GRB_DoubleAttr_Runtime);

		if(mlin->get(GRB_IntAttr_Status) != 2){
			cout << "\nError: model could not be solved. Gurobi Status " << mlin->get(GRB_IntAttr_Status) << endl;
			exit(2);
		}

		iter_num += 1;
		new_val = mlin->get(GRB_DoubleAttr_ObjVal);

		if(abs((new_val - old_val)/(old_val + eps_main)) < stall_tol)
			iter_stall ++;
		else{
			iter_stall = 0;
			if((new_val - old_val)/abs(old_val) < -stall_tol && mlin->get(GRB_IntAttr_ModelSense) > 0 ){
				cout << "\nWarning!! Numerical instability, finishing algorithm" << endl;
				break;
			}
			else if((new_val - old_val)/abs(old_val) > stall_tol && mlin->get(GRB_IntAttr_ModelSense) < 0){
				cout << "\nWarning!! Numerical instability, finishing algorithm" << endl;
				break;
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

	cout << endl << fullfilename ;
	printf(" INFO: | %2.12f | %2.12f | %6d | %6d | %6d | %6d | %5d | %3.2f | %3.2f | %1.2f | %3.2f | %3.2f | %3.2f\n\n",
					RLT_val, old_val, counts[0], counts[1], counts[2], counts[3], iter_num, run_time, gurobi_time, gurobi_time/run_time, pre_time, cut_time, post_time );

	delete[] const_number;
	delete[] fullx;
	vectoX.clear();
	A.clear();
	b.clear();
	c.clear();

	for(int i = 0; i< n; i++)
		delete[] Xtovec[i];
	
	delete[] Xtovec;
	
	delete[] counts;

	delete[] x;
	for(int i = 0; i< n; i++){
		delete[] X[i];
	}
	delete[] X;
	delete[] varlist;
	delete m;
	delete mlin;
	delete mlin_base;
	delete env;

  	return 0;
}

void computeBasis(GRBModel *m, vector<RowVectorXd> A, vector<double> b, GRBVar *x, GRBVar **X, int n, int ** Xtovec, int N, int M, MatrixXd *out_Abasic, VectorXd *out_bbasic, int *const_number){
	
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
			const_number[count] = j;
			count++;
		}
	}
		
	if(count < N){
		for(int i=0; i < n; i++){
			if(x[i].get(GRB_IntAttr_VBasis) == -1){
				Abasic(count,i) = -1;
				bbasic(count) = -x[i].get(GRB_DoubleAttr_LB);
				const_number[count] = -1; //to flag bound
				count++;
			}
			else if( x[i].get(GRB_IntAttr_VBasis) == -2){
				Abasic(count,i) = 1;
				bbasic(count) = x[i].get(GRB_DoubleAttr_UB);
				const_number[count] = -1;
				count++;
			}
				
			if(count == N) break;
		}
			
		for(int l=0; l<n; l++){
			for(int k=0; k<l+1; k++){
				if(X[k][l].get(GRB_IntAttr_VBasis) == -1){
					Abasic(count, Xtovec[k][l]) = -1;
					bbasic(count) = -X[k][l].get(GRB_DoubleAttr_LB);
					const_number[count] = -1;
					count++;
				}
				else if(X[k][l].get(GRB_IntAttr_VBasis) == -2){
					Abasic(count, Xtovec[k][l]) = 1;
					bbasic(count) = X[k][l].get(GRB_DoubleAttr_UB);
					const_number[count] = -1;
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
