#include "cutGenerators.h"

double eps = 1E-8;
double eps_check = 1E-6;
double stepback=1E-5; //what fraction to step back to avoid infeas cuts
bool doProxy = true;

Vector2d computeRoots(double a, double b, double c);
Vector3d linearFormProd(Vector2d c1, Vector2d c2);
double myMatrixInnerProduct(MatrixXd m1, MatrixXd m2, int N);

bool isInSOC(double x, double y, double z, double interiortol, double *out_violation){
	// check if ||x,y|| < z
	if (z > interiortol && x*x + y*y - z*z < -interiortol){
		*out_violation = (z*z - x*x - y*y)/max(z*z, x*x + y*y);
		return true;
	 }
	return false;
}

void generalizedminorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix,
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation){

	double PDtol=1E-8;
	double interiortol = 1E-4; //what is considered to be in the interior of the cone
	bool safe_strengthening = true;

	int counter = 0;
	vector<tuple<int,int, int, int, double, char, int, int>> minor_violation_tuples;
	double maxviol = 0;

	for(int i=0; i<n; i++)
		for(int j=i+1; j<n + 1; j++)
			for(int k=0; k<n; k++)
				for(int l=k+1; l<n+1; l++){

					double a = solX_matrix(i,k);
					double b = solX_matrix(i,l);
					double c = solX_matrix(j,k);
					double d = solX_matrix(j,l);

					double normviol;

					//Check (10a)
					double x_check = (b + c)/2;
					double y_check = (a - d)/2;
					if( i != l && j != k && isInSOC(x_check, y_check, (a + d)/2, interiortol, &normviol) ){ //neither b nor c are diagonal
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 1, 0));
					}
					else if( i == l && isInSOC(x_check, y_check, (b - c)/2, interiortol, &normviol) ){ //b is diagonal
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 0, 1));
					}
					else if( j == k && isInSOC(x_check, y_check, - (b - c)/2, interiortol, &normviol) ){ //c is diagonal
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 0, -1));
					}
					else if( i != l && j != k && i != k && j!= l) { // none of them is diagonal
						if( isInSOC(x_check, y_check, (a + d)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 1, 0));
						}
						if( isInSOC(x_check, y_check, (b - c)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 0, 1));
						}

						//new extremes
						if( isInSOC(x_check, y_check, -(a + d)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', -1, 0));
						}
						if( isInSOC(x_check, y_check, -(b - c)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a', 0, -1));
						}

					}

					//Check (10b)
					x_check = (a + d)/2;
					y_check = (b - c)/2;
					if( (i == l || j == k) && isInSOC(x_check, y_check, (b + c)/2, interiortol, &normviol) ){ //either b or c are diagonal
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 1, 0));
					}
					else if( i == k && j !=l && isInSOC(x_check, y_check, (a - d)/2, interiortol, &normviol) ){ //a is diagonal but not d
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 0, 1));
					}
					else if( i != k && j ==l && isInSOC(x_check, y_check, -(a - d)/2, interiortol, &normviol) ){ //d is diagonal but not a
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 0, -1));
					}
					else if( i != l && j != k && i != k && j!= l) { // none of them is diagonal
						if( isInSOC(x_check, y_check, (b + c)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 1, 0));
						}
						if( isInSOC(x_check, y_check, (a - d)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 0, 1));
						}

						if( isInSOC(x_check, y_check, -(b + c)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', -1, 0));
						}
						if( isInSOC(x_check, y_check, -(a - d)/2, interiortol, &normviol) ){
							counter++;
							minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b', 0, -1));
						}
					}
				}

	//Cut depth proxy
	if(doProxy){
		double avg_violation = 0;
		MatrixXd avg_Dir(n+1,n+1);
		avg_Dir.setZero();

		for(int dir_num=0; dir_num < N; dir_num++)
			avg_Dir = avg_Dir + dirs_matrix[dir_num];

		avg_Dir = avg_Dir/N;

		for(int vind=0; vind<counter; vind++){

			int i = get<0>(minor_violation_tuples[vind]);
			int j = get<1>(minor_violation_tuples[vind]);
			int k = get<2>(minor_violation_tuples[vind]);
			int l = get<3>(minor_violation_tuples[vind]);

			char set_type = get<5>(minor_violation_tuples[vind]);

			int mu1 = get<6>(minor_violation_tuples[vind]);
			int mu2 = get<7>(minor_violation_tuples[vind]);

			MatrixXd D = avg_Dir;

			Vector2d linA(D(i,k),solX_matrix(i,k)); //a = lambda*D11 + X11
			Vector2d linB(D(i,l),solX_matrix(i,l)); //b = lambda*D12 + X12
			Vector2d linC(D(j,k),solX_matrix(j,k)); //c = lambda*D21 + X21
			Vector2d linD(D(j,l),solX_matrix(j,l)); //d = lambda*D22 + X22

			double steplength = std::numeric_limits<double>::infinity();
			Vector2d linZ;
			Vector2d linX;
			Vector2d linY;

			if(set_type == 'a'){
				linZ = mu1*(linA + linD)/2 + mu2*(linB - linC)/2;
				linX = (linB + linC)/2;
				linY = (linA - linD)/2;
			}
			else if(set_type == 'b'){
				linZ = mu1*(linB + linC)/2 + mu2*(linA - linD)/2;
				linX = (linA + linD)/2;
				linY = (linB - linC)/2;
			}
			if(abs(linZ(0)) > eps && -linZ(1)/linZ(0) > eps)
				steplength = min(steplength, -linZ(1)/linZ(0)); //steplengths left-hand side == 0

			Vector3d quad = linearFormProd(linZ,linZ) - linearFormProd(linX,linX) - linearFormProd(linY,linY);
			Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

			if(roots(0) > eps)
				steplength = min(steplength, roots(0));

			if(roots(1) > eps)
				steplength = min(steplength, roots(1));

			get<4>(minor_violation_tuples[vind]) = steplength;
		}
	}

	sort(begin(minor_violation_tuples), end(minor_violation_tuples),
			  [](tuple<int, int, int, int, double, char, int, int > const &t1, tuple<int, int, int, int, double, char, int, int> const &t2)
			     {return get<4>(t1) > get<4>(t2);} );

	vector<RowVectorXd> pi_all;
	vector<double> pirhs_all;
	vector<double> violation_all;
	int minorscount = min(counter, max_cuts);

	for(int vind=0; vind < minorscount; vind++){

		int i = get<0>(minor_violation_tuples[vind]);
		int j = get<1>(minor_violation_tuples[vind]);
		int k = get<2>(minor_violation_tuples[vind]);
		int l = get<3>(minor_violation_tuples[vind]);

		char set_type = get<5>(minor_violation_tuples[vind]);

		int mu1 = get<6>(minor_violation_tuples[vind]);
		int mu2 = get<7>(minor_violation_tuples[vind]);

		//find intersection points
		RowVectorXd Beta(N);
		Beta.setZero();
		MatrixXd lamDm[N];

		int infflags[N]; // 0 = finite step, 1 = infinite step to strengthen, 2 = infinite step not to strengthen
		for(int idx=0; idx<N; idx++)
			infflags[idx] = 0;

		for( int dir_num=0; dir_num < N; dir_num++){
			MatrixXd D = dirs_matrix[dir_num];
			double a = D(i,k);
			double b = D(i,l);
			double c = D(j,k);
			double d = D(j,l);

			double cone_viol;
			double steplength;
			//check if direction is also in the cone
			if(set_type == 'a' && isInSOC((b+c)/2, (a-d)/2, mu1*(a+d)/2 + mu2*(b-c)/2, 0, &cone_viol)){
				steplength = std::numeric_limits<double>::infinity();
				//if additionally is in the *interior* of the cone, mark for strengthening
				if(isInSOC((b+c)/2, (a-d)/2, mu1*(a+d)/2 + mu2*(b-c)/2, eps, &cone_viol)){
					infflags[dir_num] = 1;
				}
				else{
					infflags[dir_num] = 2;
				}
			}
			else if(set_type == 'b' && isInSOC((a+d)/2, (b-c)/2, mu1*(b+c)/2 + mu2*(a-d)/2, 0, &cone_viol)){
				steplength = std::numeric_limits<double>::infinity();
				//if additionally is in the *interior* of the cone, mark for strengthening
				if(isInSOC((b+c)/2, (a-d)/2, mu1*(a+d)/2 + mu2*(b-c)/2, eps, &cone_viol)){
					infflags[dir_num] = 1;
				}
				else{
					infflags[dir_num] = 2;
				}
			}
			else{
				Vector2d linA(D(i,k),solX_matrix(i,k)); //a = lambda*D11 + X11
				Vector2d linB(D(i,l),solX_matrix(i,l)); //b = lambda*D12 + X12
				Vector2d linC(D(j,k),solX_matrix(j,k)); //c = lambda*D21 + X21
				Vector2d linD(D(j,l),solX_matrix(j,l)); //d = lambda*D22 + X22

				steplength = std::numeric_limits<double>::infinity();
				Vector2d linZ;
				Vector2d linX;
				Vector2d linY;

				if(set_type == 'a'){
					linZ = mu1*(linA + linD)/2 + mu2*(linB - linC)/2;
					linX = (linB + linC)/2;
					linY = (linA - linD)/2;
				}
				else if(set_type == 'b'){
					linZ = mu1*(linB + linC)/2 + mu2*(linA - linD)/2;
					linX = (linA + linD)/2;
					linY = (linB - linC)/2;
				}
				if(abs(linZ(0)) > eps && -linZ(1)/linZ(0) > eps)
					steplength = min(steplength, -linZ(1)/linZ(0)); //steplengths making left-hand side == 0

				Vector3d quad = linearFormProd(linZ,linZ) - linearFormProd(linX,linX) - linearFormProd(linY,linY);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));

				if(steplength < std::numeric_limits<double>::infinity()){
					lamDm[dir_num] = steplength*D;
				}
				else{ //in case some numerical issue generates an infinite steplength, we don't do strengthening
					infflags[dir_num] = 2;
					safe_strengthening = false; // flag this as unsafe
					//printf("Steplength unsafe\n");
				}
			}
			if(steplength == std::numeric_limits<double>::infinity()) Beta(dir_num) = 0;
			else{
				steplength *= (1-stepback);
				Beta(dir_num) = 1.0/steplength;
			}

		}
		if(Beta.squaredNorm() < eps){
			cout << "Minor cut shows infeasible problem" << endl;
			exit(1);
		}

		//Balas' Formula
		RowVectorXd pi = Beta*Abasic;
		double pirhs = Beta.dot(bbasic) - 1;
		if ( pi.dot(solX) < pirhs ){
			pi = -pi;
			pirhs = -pirhs;
		}

		double norm = pi.lpNorm<1>();
		pi = pi/norm;
		pirhs = pirhs/norm;

		double violation = pi.dot(solX) - pirhs;

		pi_all.push_back(pi);
		pirhs_all.push_back(pirhs);
		violation_all.push_back(violation);
	}
	*out_pi = pi_all;
	*out_pirhs = pirhs_all;
	*out_violation = violation_all;
}

void principalminorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix,
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts, bool strengthen, VectorXd truesol, bool checksol,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation){

	//generate elementary minor cuts
	//to ensure separation, problem MUST have nonnegative diagonals
	double PDtol=1E-8;
	double interiortol = 1E-4; //what is considered to be in the interior of the 2x2PSD cone
	bool safe_strengthening = true;

	//find elementary violations
	int counter = 0;
	vector<tuple<int,int, double>> minor_violation_tuples;

	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n+1; j++){
			if(solX_matrix(i,i) > interiortol && solX_matrix(j,j) > interiortol && solX_matrix(i,i)*solX_matrix(j,j) - solX_matrix(i,j)*solX_matrix(i,j) > interiortol){
				double normviol = (solX_matrix(i,i)*solX_matrix(j,j) - solX_matrix(i,j)*solX_matrix(i,j))/max(solX_matrix(i,i)*solX_matrix(j,j), solX_matrix(i,j)*solX_matrix(i,j));
				counter++;
				minor_violation_tuples.push_back(make_tuple(i,j, normviol));
			}
		}
	}

	sort(begin(minor_violation_tuples), end(minor_violation_tuples), [](tuple<int, int, double> const &t1, tuple<int, int, double> const &t2)
			     {return get<2>(t1) > get<2>(t2);} );

	vector<RowVectorXd> pi_all;
	vector<double> pirhs_all;
	vector<double> violation_all;

	int minorscount = min(counter, max_cuts);

	for(int vind=0; vind < minorscount; vind++){

		int ind1 = get<0>(minor_violation_tuples[vind]);
		int ind2 = get<1>(minor_violation_tuples[vind]);

		//find intersection points
		RowVectorXd Beta(N);
		Beta.setZero();
		MatrixXd lamDm[N];

		RowVectorXd pi;
		double pirhs;

		int infflags[N]; // 0 = finite step, 1 = infinite step with PD, 2 = infinite step without PSD
		for(int i=0; i<N; i++)
			infflags[i] = 0;

		for( int k=0; k < N; k++){
			MatrixXd D = dirs_matrix[k];
			double steplength;
			if (D(ind1,ind1) >= 0 && D(ind2,ind2) >= 0 && D(ind1,ind1)*D(ind2,ind2) - D(ind1,ind2)*D(ind1,ind2) >= 0 ){
				//direction submatrix is PSD
				steplength = std::numeric_limits<double>::infinity();

				if (D(ind1,ind1) >= eps && D(ind2,ind2) >= eps && D(ind1,ind1)*D(ind2,ind2) - D(ind1,ind2)*D(ind1,ind2) >= eps ){
					infflags[k] = 1; //just mark for strengthening the ones with positive DEFINITE matrix
				}
				else{
					infflags[k] = 2;
				}
			}
			else{
				Vector2d linA(D(ind1,ind1),solX_matrix(ind1,ind1)); //a = X11 + lambda*D11
				Vector2d linB(D(ind1,ind2),solX_matrix(ind1,ind2)); //b = X12 + lambda*D12
				Vector2d linC(D(ind2,ind2),solX_matrix(ind2,ind2)); //c = X22 + lambda*D22

				steplength = std::numeric_limits<double>::infinity();

				if(abs(linA(0)) > eps && -linA(1)/linA(0) > eps)
					steplength = min(steplength, -linA(1)/linA(0)); //steplengths making a = 0

				if(abs(linC(0)) > eps && -linC(1)/linC(0) > eps)
					steplength = min(steplength, -linC(1)/linC(0)); //steplengths making c = 0

				Vector3d quad = linearFormProd(linA,linC) - linearFormProd(linB,linB);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));

				if(steplength < std::numeric_limits<double>::infinity()){
					lamDm[k] = steplength*D;
				}
				else{
					infflags[k] = 2;
					safe_strengthening = false; // flag this as unsafe
				}
			}
			if(steplength == std::numeric_limits<double>::infinity()) Beta(k) = 0;
			else{
				steplength *= (1-stepback);
				Beta(k) = 1.0/steplength;
			}
		}

		if(Beta.squaredNorm() < eps){
			cout << "Minor cut shows infeasible problem" << endl;
			exit(1);
		}
		if(checksol){
			pi = Beta*Abasic;
			pirhs = Beta.dot(bbasic) - 1;
			if ( pi.dot(solX) < pirhs ){
				pi = -pi;
				pirhs = -pirhs;
			}

			double norm = pi.lpNorm<1>();
			pi = pi/norm;
			pirhs = pirhs/norm;

			if(pi.dot(truesol) - pirhs > eps_check){
				cout << "Error! Cut before strengthening cuts off given sol by " << pi.dot(truesol) - pirhs << endl;
			}
		}

		//STRENGTHENING STARTS HERE
		if(strengthen && safe_strengthening){
			for(int k=0; k<N; k++){
				if(infflags[k] == 1){
					MatrixXd D = dirs_matrix[k];
					double steplength = std::numeric_limits<double>::infinity();
					for(int m=0; m<N; m++){ //loop over all finite steps
						if(infflags[m] != 0) continue;

						double new_steplength = std::numeric_limits<double>::infinity();

						//Reminder! This steplength should be negative!
						if(abs(D(ind1,ind1)) > eps && (lamDm[m])(ind1,ind1)/D(ind1,ind1) < -eps )
							new_steplength = min(new_steplength, (lamDm[m])(ind1,ind1)/D(ind1,ind1) ); //steplengths making D11 = 0

						if(abs(D(ind2,ind2)) > eps && (lamDm[m])(ind2,ind2)/D(ind2,ind2) < -eps )
							new_steplength = min(new_steplength, (lamDm[m])(ind2,ind2)/D(ind2,ind2) ); //steplengths making D22 = 0

						double a = D(ind1,ind1)*D(ind2,ind2) - D(ind1,ind2)*D(ind1,ind2);
						double b = - (lamDm[m])(ind1,ind1)*D(ind2,ind2) - (lamDm[m])(ind2,ind2)*D(ind1,ind1) + 2*(lamDm[m])(ind1,ind2)*D(ind1,ind2);
						double c = (lamDm[m])(ind1,ind1)*(lamDm[m])(ind2,ind2) - (lamDm[m])(ind1,ind2)*(lamDm[m])(ind1,ind2);

						Vector2d roots = computeRoots(a, b, c);
						//cout << "Roots for potential strengthening " << roots.transpose() << endl;

						if(roots(0) < -eps)
							new_steplength = min(new_steplength, roots(0));

						if(roots(1) < -eps)
							new_steplength = min(new_steplength, roots(1));

						new_steplength *= (1+stepback);
						//if(new_steplength > steplength)
						//		cout << "Steplength in "<< k << " weakened due " << m << " to " << new_steplength << endl;
						steplength = min(steplength, new_steplength);
					}//for loop

					double old_beta = Beta(k);
					Beta(k) = 1/(steplength);

					if(checksol){
						pi = Beta*Abasic;
						pirhs = Beta.dot(bbasic) - 1;
						if ( pi.dot(solX) < pirhs ){
							pi = -pi;
							pirhs = -pirhs;
						}

						if(pi.dot(truesol) > pirhs){
							cout << "Error! Cut became invalid when strengthening "<< k << " by " << pi.dot(truesol) - pirhs << endl;
							cout << "Steplength given by strengthening was " << steplength <<  " New beta(k) " << Beta(k) << " Old beta(k) " << old_beta << endl;
							checksol = false;
						}
					}
				}
			}
		}
		//Balas' Formula
		pi = Beta*Abasic;
		pirhs = Beta.dot(bbasic) - 1;
		if ( pi.dot(solX) < pirhs ){
			pi = -pi;
			pirhs = -pirhs;
		}
		double norm = pi.lpNorm<1>();
		pi = pi/norm;
		pirhs = pirhs/norm;
		double violation = pi.dot(solX) - pirhs;

		pi_all.push_back(pi);
		pirhs_all.push_back(pirhs);
		violation_all.push_back(violation);
	}
	*out_pi = pi_all;
	*out_pirhs = pirhs_all;
	*out_violation = violation_all;
}

void outerpsd(VectorXd solX, MatrixXd solX_matrix, int n, int N, int **Xtovec,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation){

	double eigtol = -1E-6;
	EigenSolver<MatrixXd> es(solX_matrix);
	VectorXd eigvals = es.eigenvalues().real();
	MatrixXd V = es.eigenvectors().real();

	vector<RowVectorXd> pi_all;
	vector<double> pirhs_all;
	vector<double> violation_all;

	for(int i = 0; i < n + 1; i++){
		if(eigvals(i) < eigtol){ //eigenvalue negative enough
			VectorXd rvec = V.col(i);
			RowVectorXd pi(N);
			pi.setZero();

			for(int j=0; j < n; j++){
				pi(j) = -2*rvec(j+1)*rvec(0);

				for(int k=0; k < j+1; k++){
					if(k != j)
						pi(Xtovec[k][j]) = -2*rvec(k+1)*rvec(j+1);
					else
						pi(Xtovec[k][j]) = -rvec(k+1)*rvec(k+1);
				}
			}
			double pirhs = rvec(0)*rvec(0);

			double norm = pi.lpNorm<1>();
			pi = pi/norm;
			pirhs = pirhs/norm;

			double violation = pi.dot(solX) - pirhs;

			pi_all.push_back(pi);
			pirhs_all.push_back(pirhs);
			violation_all.push_back(violation);
		}
	}
	*out_pi = pi_all;
	*out_pirhs = pirhs_all;
	*out_violation = violation_all;
}

void shiftedconeeymcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix,
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, bool strengthen,
		RowVectorXd *out_pi, double *out_pirhs, double *out_violation){

	double bignum = 1E10;
	EigenSolver<MatrixXd> es(solX_matrix);

	VectorXd eigvals = es.eigenvalues().real();
	MatrixXd V = es.eigenvectors().real();

	vector<tuple<double, VectorXd>> eigenDec;

	for(int i = 0; i < n+1; i++)
		eigenDec.push_back(make_tuple(eigvals(i), V.col(i)));

	sort(begin(eigenDec), end(eigenDec),
			  [](tuple<double, VectorXd> const &t1, tuple<double, VectorXd> const &t2)
			     {return get<0>(t1) > get<0>(t2);} );

	MatrixXd tempmat = solX_matrix - get<0>(eigenDec[0])*get<1>(eigenDec[0])*(get<1>(eigenDec[0]).transpose());

	double radius = tempmat.norm();

	RowVectorXd pi;
	double pirhs;

	if(get<0>(eigenDec[0]) > eps){
		if(get<0>(eigenDec[1]) < -eps){ //cut is given by opf halfspace, can replace rhs with small pos num
			cout << "HALFSPACE" << endl;

			MatrixXd tempmat2 = get<0>(eigenDec[0])*get<1>(eigenDec[0])*get<1>(eigenDec[0]).transpose();
			pirhs = myMatrixInnerProduct(tempmat,tempmat2, n+1) - tempmat(0, 0); //The [0][0] part of solXmatrix is a constant

			MatrixXd temp3 = tempmat.diagonal().asDiagonal();
			MatrixXd tempmat_doubled = 2*tempmat - temp3;
			pi = buildSolFromMatrix(tempmat_doubled, Xtovec, N, n).transpose();
		}

		else{
			//cout << "CONE CUT" << endl;
			tempmat = get<0>(eigenDec[0])*get<1>(eigenDec[0])*(get<1>(eigenDec[0]).transpose());
			double ratio = get<0>(eigenDec[0])/get<0>(eigenDec[1]);

			MatrixXd C = tempmat + ratio*(solX_matrix - tempmat);

			double q = ratio*radius; //Shifted ball has centre C and radius q
			double Cfro = C.norm();
			double Csqs = Cfro*Cfro - q*q;
			double innerXC = myMatrixInnerProduct(solX_matrix,C, n+1);

			MatrixXd Z3 = Cfro*solX_matrix - innerXC*C/Cfro;

			RowVectorXd Beta(N);
			Beta.setZero();

			int infflags[N];
			for(int i=0; i<N; i++)
				infflags[i] = 0;

			bool testflag = true;
      		MatrixXd lamDm[N];

			for(int i=0; i<N; i++){
				MatrixXd D = dirs_matrix[i];

        		//check for finite step length
				double innerDC = myMatrixInnerProduct(D,C, n+1);
				double r1 = innerDC*q/(Cfro*sqrt(Csqs));
				MatrixXd auxMatrix = D - innerDC*C/(Cfro*Cfro);
				double d1 = auxMatrix.norm();
        		double steplength = 0;

				if(innerDC < eps || d1 + eps > r1){
					//finite intersection
					MatrixXd Z4 = Cfro*D - innerDC*C/Cfro;
					double quada = (q*q*innerDC*innerDC)/Csqs - Z4.squaredNorm();
					double quadb = 2*(q*q)*innerXC*innerDC/Csqs - 2*myMatrixInnerProduct(Z3,Z4,n+1);
					double quadc = (q*q*innerXC*innerXC)/Csqs - Z3.squaredNorm();

					Vector2d roots = computeRoots(quada, quadb, quadc);
					if( min(roots(0), roots(1)) > 0 )
						steplength = min(roots(0), roots(1));
					else
						steplength = max(roots(0), roots(1));

					steplength *= (1-stepback);
					MatrixXd testpt = solX_matrix + steplength*D;
					double testr = myMatrixInnerProduct(testpt, C, n+1)*q/(Cfro*sqrt(Csqs));

					MatrixXd auxMatrix = testpt - myMatrixInnerProduct(testpt,C, n+1)*C/(Cfro*Cfro);
					double testd = auxMatrix.norm();

					if((testr - testd)/testd < -eps ){
						cout << "STEPLENGTH ERROR " << (testr-testd)/testd << " REVERTING ---------" << endl;
						steplength = radius/D.norm();
					}
					lamDm[i] = steplength*D; //for strengthening
					testflag=0;
				}
				else{
					//INTERSECTION AT INFINITY
					steplength = std::numeric_limits<double>::infinity();
					infflags[i] = 1;
				}
				Beta(i) = 1/steplength;
			}

			//strengthen infinite steps
			if(strengthen){
				for(int i=0; i<N; i++){
					if(infflags[i] == 1){
						double steplength = std::numeric_limits<double>::infinity();
						MatrixXd D = dirs_matrix[i];

						for(int m=0; m<N; m++){
							if(infflags[m] != 0) continue;

							double innerLC = myMatrixInnerProduct(lamDm[m],C, n+1);
							MatrixXd Z1 = Cfro*lamDm[m] - innerLC*C/Cfro;

							double innerDC = myMatrixInnerProduct(D,C,n+1);
							MatrixXd Z2 = innerDC*C/Cfro-Cfro*D;
							double quada = (q*q*innerDC*innerDC)/Csqs - Z2.squaredNorm();
							double quadb = -2*(q*q)*innerLC*innerDC/Csqs - 2*myMatrixInnerProduct(Z1,Z2,n+1);
							double quadc = (q*q*innerLC*innerLC)/Csqs - Z1.squaredNorm();

							double new_steplength = std::numeric_limits<double>::infinity();
							Vector2d roots = computeRoots(quada, quadb, quadc);

							if(min(roots(0),roots(1)) >= 0 )
								cout << "Warning! Tilt error, positive steps" << endl ;//leave step at infinity

							if( roots(0) < -eps )
								new_steplength = min(roots(0), new_steplength);

							if( roots(1) < -eps )
								new_steplength = min(roots(1), new_steplength);

							new_steplength *= (1+stepback); //the safe stepback is to make the step more negative

							steplength = min(steplength, new_steplength); //keep the most negative among all finite steps
							}

							Beta(i) = 1/steplength;
						} // if direction has finite steps
				} //for loop in directions
			} //if strengthen
			pi = Beta*Abasic;
			pirhs = Beta.dot(bbasic) - 1;
		}
	}
	else{
		cout << "OA CUT" << endl;
		// <Xsol, X> <= 0
		pirhs = -1; //The [0][0] part of solXmatrix is 1
		MatrixXd temp3 = solX_matrix.diagonal().asDiagonal();
		MatrixXd tempmat = 2*solX_matrix - temp3;
		pi = buildSolFromMatrix(tempmat, Xtovec, N, n).transpose();
	}

	if ( pi.dot(solX) < pirhs ){
		pi = -pi;
		pirhs = -pirhs;
	}

	double norm = pi.lpNorm<1>();
	pi = pi/norm;
	pirhs = pirhs/norm;
	double violation = pi.dot(solX) - pirhs;

	*out_pi = pi;
	*out_pirhs = pirhs;
	*out_violation = violation;
}

void eymcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix,
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic,
		RowVectorXd *out_pi, double *out_pirhs, double *out_violation){

	//given polyhedron defined by Ax<=b, generate EYM
	//intersection ball cut centered at (or near) x
	//return coeffs for the cut pi(x) <= pirhs

	VectorXcd eigvalsC = solX_matrix.eigenvalues();

	if(eigvalsC.imag().squaredNorm() > eps){
		cout << "Error: Symmeytric matrix with complex eigenvalues" << endl;
		exit(1);
	}

	VectorXd eigvals = eigvalsC.real();

	sort(eigvals.data(),eigvals.data()+eigvals.size(), [](double a, double b){return a > b;} );

	double radius;
	if( eigvals[0] < 0 ){
		radius = (1-stepback)*solX_matrix.norm();
	}
	else{
		eigvals[0] = 0; //I do this to ignore the first eigenvalue
		radius = (1-stepback)*eigvals.norm();
	}

	RowVectorXd Beta(N);
	Beta.setZero();

	for(int k=0; k<N; k++){
		MatrixXd D = dirs_matrix[k];
		Beta(k)= D.norm()/radius; //Inverse step length
	}

	RowVectorXd pi = Beta*Abasic;
	double pirhs = Beta.dot(bbasic) - 1;

	if ( pi.dot(solX) < pirhs ){
		pi = -pi;
		pirhs = -pirhs;
	}

	double norm = pi.lpNorm<1>();
	pi = pi/norm;
	pirhs = pirhs/norm;

	double violation = pi.dot(solX) - pirhs;

	*out_pi = pi;
	*out_pirhs = pirhs;
	*out_violation = violation;
}


void findRays(MatrixXd Abasic, int **Xtovec, int n, int N, vector<VectorXd> *out_dirs, vector<MatrixXd> *out_dirs_matrix){
	MatrixXd inverseB = Abasic.inverse();
	vector<VectorXd> dirs(N);
	vector<MatrixXd> dirs_matrix(N);

	for(int i=0; i<N; i++){
		dirs[i] = -inverseB.col(i);
	}

	for(int k=0; k<N; k++){
		MatrixXd D(1+n, 1+n);
		D.setZero();

		VectorXd dir_vec = dirs[k];

		for(int j=0 ; j < n; j++){
			double val = dir_vec(j);
			D(j+1,0) = val;
			D(0,j+1) = val;

			for(int i=0; i<j+1; i++){
				val = dir_vec(Xtovec[i][j]);
				D(i+1,j+1) = val;
				D(j+1,i+1) = val;
			}
		}
		dirs_matrix[k] = D;
	}
	*out_dirs = dirs;
	*out_dirs_matrix = dirs_matrix;
}

MatrixXd buildMatrixSol(VectorXd solX, int **Xtovec, int n){ //We assume solX = [x, vect(X)]

	MatrixXd solX_matrix(1+n, 1+n);
	solX_matrix.setZero();

	solX_matrix(0,0) = 1;

	for(int j=0; j<n; j++){

		solX_matrix(j+1,0) = solX(j);
		solX_matrix(0,j+1) = solX(j);

		for(int i=0; i < j+1; i++){
			solX_matrix(i+1,j+1) = solX(Xtovec[i][j]);
			solX_matrix(j+1,i+1) = solX(Xtovec[i][j]);
		}
	}
	return solX_matrix;
}

VectorXd buildSolFromMatrix( MatrixXd solX_matrix, int **Xtovec, int N, int n){
	VectorXd solX(N);
	solX.setZero();

	for(int j=0; j < n; j++){
		solX(j) = solX_matrix(j+1,0);

		for(int i=0; i<j+1; i++){
			solX(Xtovec[i][j]) = solX_matrix(i+1,j+1);
		}
	}
	return solX;
}

Vector3d linearFormProd( Vector2d c1, Vector2d c2){
	double a = c1(0)*c2(0);
	double b = c1(0)*c2(1) + c1(1)*c2(0);
	double c = c1(1)*c2(1);

	Vector3d result(a,b,c);
	return result;
}

Vector2d computeRoots(double a, double b, double c){
	//compute real roots of ax^2 +bx + c
	if(abs(a) < eps){ //quad formula is approx linear, use more stable lin eq.
		if(abs(b) < eps){
			Vector2d result(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
			return result;
		}
		else{
			Vector2d result(c/(-b), c/(-b));
			return result;
		}
	}
	double disc = b*b - 4*a*c;
	double disc_test = disc/max(b*b,abs(4*a*c));

	if(disc_test < -eps){
		Vector2d result(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
		return result;
	}
	else if(abs(disc_test) < eps){ //discriminant nearly 0 (scaled)
		Vector2d result(-b/(2*a), -b/(2*a));
		return result;
	}
	else{
		double disc_sqrt = sqrt(disc);
		Vector2d result((-b+disc_sqrt)/(2*a) , (-b-disc_sqrt)/(2*a));
		return result;
	}
}
int checkifparallel(RowVectorXd pi, double pirhs, vector<RowVectorXd> A, vector<double> b, int M, vector<double> rowNorms){

	double normpi = pi.norm();

	for(int i=0; i < M; i++){
		RowVectorXd v = A[i];
		double normv = rowNorms[i];
		double inner = pi.dot(v);
		if(normpi*normv - abs(inner) < eps && inner > 0 and b[i]/normv <= pirhs/normpi)
			return 1;
	}
	return 0;
}
double myMatrixInnerProduct(MatrixXd m1, MatrixXd m2, int N){
	double prod = 0;
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			prod += m1(i,j)*m2(i,j);

	return prod;
}
