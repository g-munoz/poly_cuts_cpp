#include "formulationHandler.h"

void boundTightening(GRBModel *m, GRBVar* varlist, int n, map<string,int> varNameToNumber, map<int,string> varNumberToName){

	bool boundImproved = false;
	GRBVar *x;
	GRBVar **X;
	GRBModel *bound_model = linearize(m, varNameToNumber, varNumberToName, false, &x, &X);
		
	for(int i=0; i < n; i++){
		GRBLinExpr obj = GRBLinExpr(x[i]);

		bound_model->setObjective(obj, GRB_MAXIMIZE);
		bound_model->update();
		bound_model->getEnv().set(GRB_IntParam_OutputFlag,0);
		bound_model->optimize();
			
		if(bound_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			double newUB = bound_model->get(GRB_DoubleAttr_ObjVal);
			if(newUB < x[i].get(GRB_DoubleAttr_UB)){
				boundImproved = true;
				cout << "Upper Bound for " << x[i].get(GRB_StringAttr_VarName) << " changed from " << x[i].get(GRB_DoubleAttr_UB) << " to " << newUB << endl;
				varlist[i].set(GRB_DoubleAttr_UB, newUB);
			}
		}
				
		bound_model->setObjective(obj, GRB_MINIMIZE);
		bound_model->update();
		bound_model->optimize();
			
		if(bound_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			double newLB = bound_model->get(GRB_DoubleAttr_ObjVal);
			if(newLB > x[i].get(GRB_DoubleAttr_LB)){
				boundImproved = true;
				cout << "Lower Bound for " << x[i].get(GRB_StringAttr_VarName) << " changed from " << x[i].get(GRB_DoubleAttr_LB) << " to " << newLB << endl;
				varlist[i].set(GRB_DoubleAttr_LB, newLB);
			}
		}
	}
	if(!boundImproved)
		cout << "No bound improved" << endl;
	
	m->update();
	delete bound_model;
	delete[] x;
	for(int i = 0; i< n; i++){
		delete[] X[i];
	}
	delete[] X;
}


GRBModel* linearize(GRBModel *m, map<string,int> varNameToNumber, map<int,string> varNumberToName, bool wRLT,
	GRBVar **out_x, GRBVar ***out_X){
	// Linearize all quadratic monomials and return new model
	// We assume all variables are named xi
	
	// wRLT determines if weak RLT is used or not.
	
	GRBModel *mlin = new GRBModel(*m);
	GRBVar* x = mlin->getVars(); // save the original variables
	int numQConstrs = mlin->get(GRB_IntAttr_NumQConstrs);
	int n = mlin->get(GRB_IntAttr_NumVars);
	
	GRBVar **X = new GRBVar*[n]; //save the linearized monomials
	for(int i = 0; i< n; i++){
		X[i] = new GRBVar[n];
	}

	int isInitialized[n][n];
	for(int i = 0; i< n; i++)
		for(int j = 0; j< n; j++)
			isInitialized[i][j] = 0;

	GRBQuadExpr obj = mlin->getObjective();
	GRBLinExpr linear_obj = obj.getLinExpr();

	double coeffs[obj.size()];
	GRBVar vars[obj.size()];
	char name[50];

	for(int j=0; j < obj.size(); j++){
			
		string name1 = obj.getVar1(j).get(GRB_StringAttr_VarName);
		string name2 = obj.getVar2(j).get(GRB_StringAttr_VarName);
			
		int var1 = varNameToNumber[name1];
		int var2 = varNameToNumber[name2];

		if(var1 == var2){
			if(isInitialized[var1][var2] == 0){ //if not added already
				sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
				X[var1][var2] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var1][var2] = 1;
				mlin->update();
			}
			coeffs[j] = obj.getCoeff(j);
			vars[j] = X[var1][var2];
		}

		else if(var1 < var2){
			if(isInitialized[var1][var2] == 0){
				sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
				X[var1][var2] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var1][var2] = 1;
				mlin->update();
			}
			coeffs[j] = obj.getCoeff(j);
			vars[j] = X[var1][var2];
		}

		else{
			if(isInitialized[var2][var1] == 0){
				sprintf(name, "X(%s,%s)", name2.c_str(), name1.c_str());
				X[var2][var1] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var2][var1] = 1;
				mlin->update();
			}
			coeffs[j] = obj.getCoeff(j);
			vars[j] = X[var2][var1];
		}
	}
	
	mlin->update();
	linear_obj.addTerms(coeffs, vars, obj.size());
	mlin->setObjective(linear_obj);
	
	for(int i=0; i < numQConstrs; i++){ // linearize quadratic constraints
		
		GRBQConstr q = mlin->getQConstrs()[0];	
		GRBQuadExpr quad = mlin->getQCRow(q);
		GRBLinExpr linear = quad.getLinExpr();

		double coeffsQ[quad.size()];
		GRBVar varsQ[quad.size()];
		
		for(int j=0 ; j< quad.size(); j++){
			
			string name1 = quad.getVar1(j).get(GRB_StringAttr_VarName);
			string name2 = quad.getVar2(j).get(GRB_StringAttr_VarName);
			
			int var1 = varNameToNumber[name1];
			int var2 = varNameToNumber[name2];

			if(var1 == var2){
				if(isInitialized[var1][var2] == 0){ //if not added already
					sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
					X[var1][var2] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var1][var2] = 1;
					mlin->update();
				}

				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = X[var1][var2];
			}

			else if(var1 < var2){
				if(isInitialized[var1][var2] == 0){
					sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
					X[var1][var2] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var1][var2] = 1;
					mlin->update();
				}
				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = X[var1][var2];
			}

			else{
				if(isInitialized[var2][var1] == 0){
					sprintf(name, "X(%s,%s)", name2.c_str(), name1.c_str());
					X[var2][var1] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var2][var1] = 1;
					mlin->update();
				}
				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = X[var2][var1];
			}
		}
		
		linear.addTerms(coeffsQ, varsQ, quad.size());
		mlin->update() ;

		string name_str = q.get(GRB_StringAttr_QCName) + "_lin" ;
		mlin->addConstr(linear , q.get(GRB_CharAttr_QCSense), q.get(GRB_DoubleAttr_QCRHS), name_str) ;
		mlin->remove(mlin->getQConstrs()[0]);
		mlin->update();
	}

	
	for(int j = 0; j < n ; j++)
		for(int i = 0; i < j+1; i++){
			if(isInitialized[i][j] == 0){

				string namei = varNumberToName[i];
				string namej = varNumberToName[j];

				sprintf(name, "X(%s,%s)", namei.c_str(), namej.c_str());
				if(i==j)
					X[i][j] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				else
					X[i][j] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
			}
		}

	mlin->update();
	addRLTconstraints(mlin, x, X, n, wRLT);
	mlin->update();
	*out_x = x;
	*out_X = X;

	return mlin;
}

void addRLTconstraints(GRBModel *m, GRBVar* x, GRBVar** X, int n, bool wRLT){
	if(wRLT)
	    cout << "Relaxing using wRLT" << endl;
    else
        cout << "Relaxing using RLT" << endl;

	for(int j=0; j < n ; j++){
		double ubj = x[j].get(GRB_DoubleAttr_UB);
		double lbj = x[j].get(GRB_DoubleAttr_LB);
		for(int i=0; i < j+1 ; i ++){
			double ubi = x[i].get(GRB_DoubleAttr_UB);
			double lbi = x[i].get(GRB_DoubleAttr_LB);

			char name[50];
			//this is for wRLT
			if(wRLT && i == j){
				sprintf(name, "RLT4(%d,%d)", i,j);
				m->addConstr( X[i][j] <= abs(ubi*lbi), name);
			}
				
			else{
				if (lbj != - GRB_INFINITY && lbi != - GRB_INFINITY){
					sprintf(name, "RLT1(%d,%d)", i,j);
					m->addConstr( X[i][j] - lbj*x[i] - lbi*x[j] + lbi*lbj >= 0, name);
				}
				if (ubj != GRB_INFINITY && ubi != GRB_INFINITY){
					sprintf(name, "RLT2(%d,%d)", i,j);
					m->addConstr( X[i][j] - ubj*x[i] - ubi*x[j] + ubi*ubj >= 0, name);
				}
				if (lbj != - GRB_INFINITY && ubi != GRB_INFINITY){
					sprintf(name, "RLT3(%d,%d)", i,j);
					m->addConstr( X[i][j] - lbj*x[i] - ubi*x[j] + ubi*lbj <= 0, name);
				}
				if (ubj != GRB_INFINITY && lbi != - GRB_INFINITY){
					sprintf(name, "RLT4(%d,%d)", i,j);
					m->addConstr( X[i][j] - ubj*x[i] - lbi*x[j] + lbi*ubj <= 0, name);
				}
			}
		}
	}
}



void createMap(GRBVar *x, GRBVar **X, int ***out_Xtovec, vector< array<int, 2> > *out_vectoX, int n){
	// Create map between X and its vector form, to go back and forth.
	int **Xtovec = new int*[n];
	for(int i = 0; i< n; i++){
		Xtovec[i] = new int[n];
	}
	
	vector< array<int, 2> > vectoX;
	int count = 0;

	for(int j=0; j<n; j++){
		for(int i=0; i<j+1; i++){
			Xtovec[i][j] = n + count;
			array<int, 2> pair = {i,j};
			vectoX.push_back(pair);
			count += 1;
		}
	}

	*out_Xtovec = Xtovec;
	*out_vectoX = vectoX;

	return;
}


void buildAb(GRBModel *m, GRBVar *x, GRBVar **X, int **Xtovec, int n, vector<RowVectorXd> *out_A, vector<double> *out_b, vector<double> *out_c){

	// build the matrix of constraints A*y <= b
	// where y = vec([x X])

	int N = m->get(GRB_IntAttr_NumVars);
	int M = m->get(GRB_IntAttr_NumConstrs);

	vector<RowVectorXd> A;
	vector<double> b(M,0);
	vector<double> c(N,0);

	GRBConstr *constrs = m->getConstrs();
	
	// Now we go over the constraints and save the coefficients in A
	for(int j=0; j < M; j++){
		GRBLinExpr row = m->getRow(constrs[j]);
		char sense = constrs[j].get(GRB_CharAttr_Sense);
		int flip = 1;
	
		RowVectorXd row_vec(N);
		row_vec.setZero();

		if(sense == '>'){
			b[j] = -constrs[j].get(GRB_DoubleAttr_RHS);
			flip = -1;
		}
		else{
			b[j] = constrs[j].get(GRB_DoubleAttr_RHS);
		}

		for( int i = 0; i < row.size(); i++){
			double coeff = row.getCoeff(i);
			GRBVar var = row.getVar(i);
			bool found = false;

			if(var.get(GRB_StringAttr_VarName).at(0) != 'X'){
				for(int varindex = 0; varindex < n; varindex++){
					if(var.sameAs(x[varindex])){
						row_vec(varindex) = flip*coeff;
						found = true;
						break;
					}
				}
			}
			else{
				for(int l=0; l < n; l++ ){
					for(int k=0; k < l+1; k++){
						if(var.sameAs(X[k][l])){
							row_vec( Xtovec[k][l] ) = flip*coeff;
							found = true;
							break;
						}
						if(found) break;
					}
				}
			}
		}
		A.push_back(row_vec);
	}

	for(int varindex=0; varindex < n; varindex++){
		c[varindex] = x[varindex].get(GRB_DoubleAttr_Obj);
	}
	for(int l=0; l < n; l++){
		for(int k=0; k<l+1; k++){
			c[ Xtovec[k][l] ] = X[k][l].get(GRB_DoubleAttr_Obj);
		}
	}

	*out_A = A;
	*out_b = b;
	*out_c = c;

	delete[] constrs;

	return;
	
}
