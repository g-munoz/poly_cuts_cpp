#include "formulationHandler.h"


int readSol(string solfilename, int **Xtovec, map<string,int> varNameToNumber, int n, int N, VectorXd *out_vec){
	cout << "Reading sol from file "<< solfilename << endl;
	VectorXd truesol(N);
	truesol.setZero();
	ifstream infile(solfilename);
	string line;
	while (getline(infile, line))
	{
	    istringstream iss(line);
			vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

			if(results[0].compare("modelstatus") == 0 || results[0].compare("obj") == 0 || results[0].compare("nonopt") == 0 ) continue;
			if(results[0].compare("infeas") == 0 ){
				if( atof(results[2].c_str()) > 0 ){
					cout << "Error! Solution marked as infeasible" << endl;
					return 1;
				}
				continue;
			}
			auto mapSearch = varNameToNumber.find(results[0]);
			if(mapSearch == varNameToNumber.end()){
				cout << "Error! Variable " << results[0] << " not found in map" << endl;
				return 1;
			}

			int varnumber = mapSearch->second;
			truesol(varnumber) = atof(results[2].c_str());
			cout << "Variable " << results[0] << " of input sol is " << results[2] << endl;
	}
	for(int l=0; l < n; l++)
		for(int k=0; k < l+1; k++){
			if(Xtovec[k][l] == -1) continue;
			truesol(Xtovec[k][l]) = truesol(k)*truesol(l);
		}

	*out_vec = truesol;
	return 0;
}

void flushConstraints(GRBModel *mlin, int *M_ptr, int M_base, int *total_cuts_ptr, int N, MatrixXd Abasic, GRBVar *fullx, VectorXd bbasic, tuple<int,string> *const_number){

	GRBConstr *constrs = mlin->getConstrs();
	for(int l = M_base; l < (*M_ptr); l++)
		mlin->remove(constrs[l]);

	mlin->update();

	int M = M_base;
	int total_cuts = 0;

	for(int l=0; l < N; l++)
		if(get<0>(const_number[l]) >= M_base){ //if the constraint is not in the base linearization
			GRBLinExpr basis_row;
			for(int idx = 0; idx < N ; idx++)
				basis_row += Abasic(l,idx)*(fullx[idx]);
			basis_row -= bbasic[l];
			mlin->addConstr(basis_row <= 0, get<1>(const_number[l]));
			total_cuts++;
			M++;
		}

	(*total_cuts_ptr) = total_cuts;
	(*M_ptr) = M;
	delete[] constrs;
}

void boundTightening(GRBModel *m, GRBVar* varlist, int n, map<string,int> varNameToNumber, map<int,string> varNumberToName, int addAll){

	bool boundImproved = false;
	GRBVar **x;
	GRBVar ***X;
	bool **isInOriginalModel;
	GRBModel *bound_model = linearize(m, varNameToNumber, varNumberToName, false, addAll, &x, &X, &isInOriginalModel);

	for(int i=0; i < n; i++){
		GRBLinExpr obj = GRBLinExpr(*(x[i]));

		bound_model->setObjective(obj, GRB_MAXIMIZE);
		bound_model->update();
		bound_model->getEnv().set(GRB_IntParam_OutputFlag,0);
		bound_model->optimize();

		if(bound_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			double newUB = bound_model->get(GRB_DoubleAttr_ObjVal);
			if(newUB < x[i]->get(GRB_DoubleAttr_UB)){
				boundImproved = true;
				cout << "Upper Bound for " << x[i]->get(GRB_StringAttr_VarName) << " changed from " << x[i]->get(GRB_DoubleAttr_UB) << " to " << newUB << endl;
				varlist[i].set(GRB_DoubleAttr_UB, newUB);
			}
		}

		bound_model->setObjective(obj, GRB_MINIMIZE);
		bound_model->update();
		bound_model->optimize();

		if(bound_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){
			double newLB = bound_model->get(GRB_DoubleAttr_ObjVal);
			if(newLB > x[i]->get(GRB_DoubleAttr_LB)){
				boundImproved = true;
				cout << "Lower Bound for " << x[i]->get(GRB_StringAttr_VarName) << " changed from " << x[i]->get(GRB_DoubleAttr_LB) << " to " << newLB << endl;
				varlist[i].set(GRB_DoubleAttr_LB, newLB);
			}
		}
	}
	if(!boundImproved)
		cout << "No bound improved" << endl;

	m->update();
	for(int i = 0; i< n; i++){
		delete[] X[i];
		delete[] isInOriginalModel[i];
	}
	delete[] x;
	delete[] X;
	delete[] isInOriginalModel;
	delete bound_model;
}

GRBModel* linearize(GRBModel *m, map<string,int> varNameToNumber, map<int,string> varNumberToName, bool wRLT, int addAll,
	GRBVar ***out_x, GRBVar ****out_X, bool ***out_isInOriginalModel){
	// Linearize all quadratic monomials and return new model
	// We assume all variables are named xi
	// wRLT determines if weak RLT is used or not.

	GRBModel *mlin = new GRBModel(*m);
	GRBVar* x_temp = mlin->getVars(); // save the original variables
	int n = mlin->get(GRB_IntAttr_NumVars);
	GRBVar **x =  new GRBVar*[n];
	for(int i=0; i<n; i++){
		x[i] = &x_temp[i];
	}

	int numQConstrs = mlin->get(GRB_IntAttr_NumQConstrs);

	GRBVar ***X = new GRBVar**[n]; //save the linearized monomials
	for(int i = 0; i< n; i++){
		X[i] = new GRBVar*[n];
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

	double totalsum = 0.0;

	for(int j=0; j < obj.size(); j++){

		//cout << "X0,0 = " << X[0][0] << " X0,1 = " << X[0][1] << " X1,1 = " << X[1][1] << endl;

		string name1 = obj.getVar1(j).get(GRB_StringAttr_VarName);
		string name2 = obj.getVar2(j).get(GRB_StringAttr_VarName);

		int var1 = varNameToNumber[name1];
		int var2 = varNameToNumber[name2];

		coeffs[j] = obj.getCoeff(j);
		totalsum += coeffs[j];

		if(var1 == var2){
			if(isInitialized[var1][var2] == 0){ //if not added already
				sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
				X[var1][var2] = new GRBVar(); //this should be fixed
				*X[var1][var2] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var1][var2] = 1;
				mlin->update();
			}
			vars[j] = *(X[var1][var2]);
		}

		else if(var1 < var2){
			if(isInitialized[var1][var2] == 0){
				sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
				X[var1][var2] = new GRBVar();
				*X[var1][var2] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var1][var2] = 1;
				mlin->update();
			}
			//coeffs[j] = obj.getCoeff(j);
			vars[j] = *(X[var1][var2]);
		}

		else{
			if(isInitialized[var2][var1] == 0){
				sprintf(name, "X(%s,%s)", name2.c_str(), name1.c_str());
				X[var2][var1] = new GRBVar();
				*X[var2][var1] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				isInitialized[var2][var1] = 1;
				mlin->update();
			}
			//coeffs[j] = obj.getCoeff(j);
			vars[j] = *(X[var2][var1]);
		}
	}

	mlin->update();
	linear_obj.addTerms(coeffs, vars, obj.size());
	mlin->setObjective(linear_obj);

	cout << "Total sum of coeffs " << totalsum << endl;
	//exit(1);

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
					X[var1][var2] = new GRBVar();
					*X[var1][var2] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var1][var2] = 1;
					mlin->update();
				}

				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = *(X[var1][var2]);
			}

			else if(var1 < var2){
				if(isInitialized[var1][var2] == 0){
					sprintf(name, "X(%s,%s)", name1.c_str(), name2.c_str());
					X[var1][var2] = new GRBVar();
					*X[var1][var2] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var1][var2] = 1;
					mlin->update();
				}
				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = *(X[var1][var2]);
			}

			else{
				if(isInitialized[var2][var1] == 0){
					sprintf(name, "X(%s,%s)", name2.c_str(), name1.c_str());
					X[var2][var1] = new GRBVar();
					*X[var2][var1] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
					isInitialized[var2][var1] = 1;
					mlin->update();
				}
				coeffsQ[j] = quad.getCoeff(j);
				varsQ[j] = *(X[var2][var1]);
			}
		}

		linear.addTerms(coeffsQ, varsQ, quad.size());
		mlin->update() ;

		string name_str = q.get(GRB_StringAttr_QCName) + "_lin" ;
		mlin->addConstr(linear , q.get(GRB_CharAttr_QCSense), q.get(GRB_DoubleAttr_QCRHS), name_str) ;
		mlin->remove(mlin->getQConstrs()[0]);
		mlin->update();
	}

	bool **isInOriginalModel = new bool*[n];
	for(int i = 0; i < n; i++)
		isInOriginalModel[i] = new bool[n];

	for(int j = 0; j < n ; j++)
		for(int i = 0; i < j+1; i++){
			if(isInitialized[i][j] == 1){
				isInOriginalModel[i][j] = true;
				isInOriginalModel[j][i] = true;
			}
			else{
				isInOriginalModel[i][j] = false;
				isInOriginalModel[j][i] = false;
			}
		}

	// Here we add the missing quadratic terms if needed
	for(int j = 0; j < n ; j++)
		for(int i = 0; i < j+1; i++){
			if(isInitialized[i][j] == 0 && (addAll == 2 ||  (addAll == 1 && i==j) ) ){ //add diagonals if addall is 1
				string namei = varNumberToName[i];
				string namej = varNumberToName[j];
				sprintf(name, "X(%s,%s)", namei.c_str(), namej.c_str());
				if(i==j){
					X[i][j] = new GRBVar();
					*X[i][j] = mlin->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				}
				else{
					X[i][j] = new GRBVar();
					*X[i][j] = mlin->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
				}
			}
			else if(isInitialized[i][j] == 0){
				X[i][j] = nullptr;
			}
		}

	mlin->update();
	addRLTconstraints(mlin, x, X, n, wRLT);
	mlin->update();

	*out_x = x;
	*out_X = X;
	*out_isInOriginalModel = isInOriginalModel;

	return mlin;
}


GRBModel* unlinearize(GRBModel *m, GRBVar **x, GRBVar ***X, int n, int M){
	//return to non-linear model

	GRBModel *mNL = new GRBModel(m->getEnv());
	GRBVar *newX = new GRBVar[n];

	for(int i=0; i<n; i++)
		newX[i] = mNL->addVar(x[i]->get(GRB_DoubleAttr_LB), x[i]->get(GRB_DoubleAttr_UB), 0, GRB_CONTINUOUS, x[i]->get(GRB_StringAttr_VarName));

	mNL->update();

	//LinExpr obj = m->getObjective();

	GRBQuadExpr objNL = GRBQuadExpr(0);
	for(int i=0; i < n; i++){
		double coeff = x[i]->get(GRB_DoubleAttr_Obj);
		objNL.addTerm(coeff, newX[i]);

		for (int j=0; j < i+1; j++){
			coeff = 0;
			if(X[j][i] != nullptr){
				coeff = X[j][i]->get(GRB_DoubleAttr_Obj);
				objNL.addTerm(coeff, newX[j], newX[i]);
			}
		}
	}

	mNL->setObjective(objNL, m->get(GRB_IntAttr_ModelSense));
	mNL->update();

	GRBConstr *constrs = m->getConstrs();

	int M_new = m->get(GRB_IntAttr_NumConstrs);

	for(int j=0; j < M_new; j++){
		char sense = constrs[j].get(GRB_CharAttr_Sense);
		double rhs = constrs[j].get(GRB_DoubleAttr_RHS);
		string constrName = constrs[j].get(GRB_StringAttr_ConstrName);

		if (constrName.find("RLT") != std::string::npos) continue;

		bool isQuad = false;
		GRBQuadExpr rowNL = GRBQuadExpr(0);
		for(int i=0; i < n; i++){
			double coeff = m->getCoeff(constrs[j], *(x[i]));
			rowNL.addTerm(coeff, newX[i]);

			for (int k=0; k < i+1; k++){
				coeff = 0;
				if( X[k][i] != nullptr){
					coeff = m->getCoeff(constrs[j], *(X[k][i]));
					rowNL.addTerm(coeff, newX[k], newX[i]);
					if(!isQuad && abs(coeff) > 0) isQuad = true;
				}
			}
		}

		if(isQuad) mNL->addQConstr(rowNL, sense, rhs, constrName);
		else mNL->addConstr(rowNL.getLinExpr(), sense, rhs, constrName);
	}
	mNL->update();
	return mNL;
}

void projectDown(GRBModel *m, GRBVar **x, GRBVar ***X, int n, int M, bool **isInOriginalModel, bool keepRLT){

	GRBConstr *constrs = m->getConstrs();
	int *flags = new int[M]; // 1.- Leave, 2.- Remove (for RLT), 3.- Project
	GRBLinExpr *new_Constrs = new GRBLinExpr[M];
	GRBQuadExpr originalObj = m->getObjective();
	int originalModelSense = m->get(GRB_IntAttr_ModelSense);
	m->optimize();
	double objVal = m->get(GRB_DoubleAttr_ObjVal);

	for(int j=0; j < M; j++){
		GRBLinExpr row = m->getRow(constrs[j]);
		char sense = constrs[j].get(GRB_CharAttr_Sense);
		double rhs = constrs[j].get(GRB_DoubleAttr_RHS);
		string constrName = constrs[j].get(GRB_StringAttr_ConstrName);

		GRBLinExpr toProject = GRBLinExpr(0);
		GRBLinExpr toLeave = GRBLinExpr(0);
		int terms_found = 0;

		for( int i = 0; i < row.size(); i++){
			double coeff = row.getCoeff(i);
			GRBVar var = row.getVar(i);
			bool found = false;

			if(var.get(GRB_StringAttr_VarName).at(0) != 'X'){
				toLeave += coeff*var;
				continue;
			}
			else{
				for(int l=0; l < n; l++ ){
					for(int k=0; k < l+1; k++){
						if(X[k][l] != nullptr && var.sameAs(*(X[k][l]))){
							if(!isInOriginalModel[k][l]){
								toProject += coeff*var;
								terms_found++;
							}
							else{
								toLeave += coeff*var;
							}
							found = true;
							break;
						}
						if(found) break;
					}
				}
			}
		}

		new_Constrs[j] = toLeave;

		if(terms_found == 0)
			flags[j] = 1;
		else if(constrName.find("RLT") != std::string::npos){	//if a term was found, but in RLT
			if(keepRLT)
				flags[j] = 1;
			else
				flags[j] = 2;
		}
		else{
			flags[j] = 3;
			if(sense == '>')
				m->setObjective(toProject, GRB_MAXIMIZE);
			else
				m->setObjective(toProject, GRB_MINIMIZE);

			m->optimize();

			if(m->get(GRB_IntAttr_Status) != 2){
				cout << "Error, projection model could not be solved. Gurobi Status " << m->get(GRB_IntAttr_Status) << endl;
				flags[j] = 2;
			}
			else{
				new_Constrs[j] = toLeave + m->get(GRB_DoubleAttr_ObjVal) - rhs;
			}
		}
	}


	for(int j=0; j < M; j++){
		string constrName = constrs[j].get(GRB_StringAttr_ConstrName);
		double rhs = constrs[j].get(GRB_DoubleAttr_RHS);
		if(flags[j] == 1){
			continue;
		}
		else if(flags[j] == 2){
			m->remove(constrs[j]);
		}
		else{
			char sense = constrs[j].get(GRB_CharAttr_Sense);
			m->remove(constrs[j]);
			if(sense == '>')
				m->addConstr(new_Constrs[j] >= 0, constrName+"_prj");
			else
				m->addConstr(new_Constrs[j] <= 0, constrName+"_prj");
		}
	}

	if(!keepRLT){
		if(originalModelSense == GRB_MAXIMIZE){
			m->addConstr(originalObj <= objVal, "ObjCut");
		}
		else{
			m->addConstr(originalObj >= objVal, "ObjCut");
		}

		for(int l=0; l < n; l++ )
			for(int k=0; k < l+1; k++)
				if(!isInOriginalModel[k][l]){
					m->remove(*(X[k][l]));
				}
	}

	m->setObjective(originalObj, originalModelSense);
	m->update();
	return;
}

void addRLTconstraints(GRBModel *m, GRBVar** x, GRBVar*** X, int n, bool wRLT){
	if(wRLT)
		cout << "Relaxing using wRLT" << endl;
    else
    	cout << "Relaxing using RLT" << endl;

	for(int j = 0; j < n ; j++){
		double ubj = x[j]->get(GRB_DoubleAttr_UB);
		double lbj = x[j]->get(GRB_DoubleAttr_LB);
		for(int i=0; i < j+1 ; i ++){
			double ubi = x[i]->get(GRB_DoubleAttr_UB);
			double lbi = x[i]->get(GRB_DoubleAttr_LB);

			if( X[i][j] == nullptr) continue;

			char name[50];
			//this is for wRLT
			if(wRLT && i == j){
				sprintf(name, "RLT4(%d,%d)", i,j);
				m->addConstr( *(X[i][j]) <= abs(ubi*ubi), name);
			}
			else{
				if (lbj != - GRB_INFINITY && lbi != - GRB_INFINITY){
					sprintf(name, "RLT1(%d,%d)", i,j);
					m->addConstr( *(X[i][j]) - lbj*(*(x[i])) - lbi*(*(x[j])) + lbi*lbj >= 0, name);
				}
				if (ubj != GRB_INFINITY && ubi != GRB_INFINITY){
					sprintf(name, "RLT2(%d,%d)", i,j);
					m->addConstr( *(X[i][j]) - ubj*(*(x[i])) - ubi*(*(x[j])) + ubi*ubj >= 0, name);
				}
				if (lbj != - GRB_INFINITY && ubi != GRB_INFINITY){
					sprintf(name, "RLT3(%d,%d)", i,j);
					m->addConstr( *(X[i][j]) - lbj*(*(x[i])) - ubi*(*(x[j])) + ubi*lbj <= 0, name);
				}
				if (ubj != GRB_INFINITY && lbi != - GRB_INFINITY){
					sprintf(name, "RLT4(%d,%d)", i,j);
					m->addConstr( *(X[i][j]) - ubj*(*(x[i])) - lbi*(*(x[j])) + lbi*ubj <= 0, name);
				}
			}
		}
	}
}

void createMap(GRBVar **x, GRBVar ***X, int ***out_Xtovec, vector< array<int, 2> > *out_vectoX, int n){
	// Create map between X and its vector form, to go back and forth.
	int **Xtovec = new int*[n];
	for(int i = 0; i< n; i++){
		Xtovec[i] = new int[n];
		for(int j=0; j < n; j++)
			Xtovec[i][j] = -1;
	}

	vector< array<int, 2> > vectoX;
	int count = 0;

	for(int j=0; j<n; j++){
		for(int i=0; i<j+1; i++){
			if( X[i][j] == nullptr) continue;

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


void buildAb(GRBModel *m, GRBVar **x, GRBVar ***X, int **Xtovec, int n, vector<RowVectorXd> *out_A, vector<double> *out_b, RowVectorXd *out_c){

	// build the matrix of constraints A*y <= b
	// where y = vec([x X])

	int N = m->get(GRB_IntAttr_NumVars);
	int M = m->get(GRB_IntAttr_NumConstrs);

	vector<RowVectorXd> A;
	vector<double> b(M,0);
	RowVectorXd c(N);

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
					if(x[varindex] != nullptr && var.sameAs(*(x[varindex]))){
						row_vec(varindex) = flip*coeff;
						found = true;
						break;
					}
				}
			}
			else{
				for(int l=0; l < n; l++ ){
					for(int k=0; k < l+1; k++){
						if(X[k][l] != nullptr && var.sameAs(*(X[k][l]))){
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

	c.setZero();
	for(int varindex=0; varindex < n; varindex++){
		c(varindex) = x[varindex]->get(GRB_DoubleAttr_Obj);
	}
	for(int l=0; l < n; l++){
		for(int k=0; k<l+1; k++){
			if(X[k][l] == nullptr) continue;
			c( Xtovec[k][l] ) = X[k][l]->get(GRB_DoubleAttr_Obj);
		}
	}

	*out_A = A;
	*out_b = b;
	*out_c = c;

	delete[] constrs;

	return;

}
