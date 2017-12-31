#include "cutGenerators.h"

double eps = 1E-10;

bool strengthen = false;
bool doProxy = true;

Vector2d computeRoots(double a, double b, double c);
Vector3d linearFormProd( Vector2d c1, Vector2d c2);
double myMatrixInnerProduct(MatrixXd m1, MatrixXd m2, int N);


void generalizedminorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation){

	double stepback=1E-8; //what fraction to step back to avoid infeas cuts
	double PDtol=1E-8;
	double interiortol = 1E-5; //what is considered to be in the interior of the 2x2PSD cone

	//find elementary violations
	int counter = 0;
	
	vector<tuple<int,int, int, int, double, char>> minor_violation_tuples;

	double maxviol = 0;

	/*
	#test for ex2_1_1

	#print n, len(solX)
	#print Xtovec
	#optsolX = [0 for i in xrange(len(solX))]
	#optsolX[0] = 1
	#optsolX[1] = 1
	#optsolX[2] = 0
	#optsolX[3] = 1
	#optsolX[4] = 0
	#for j in xrange(n):
	#	for i in xrange(j+1):
	#		optsolX[Xtovec[i][j]] = optsolX[i]*optsolX[j]

	#print optsolX

	#raw_input()
	################################
	*/

	for(int i=0; i<n; i++)
		for(int j=i+1; j<n + 1; j++)
			for(int k=0; k<n; k++)
				for(int l=k+1; l<n+1; l++){

					double a = solX_matrix(i,k);
					double b = solX_matrix(i,l);
					double c = solX_matrix(j,k);
					double d = solX_matrix(j,l);
					
					double normviol;

					if(a > interiortol && d > interiortol && a*d - ((b + c)*(b + c)/4) > interiortol){ //this covers case (1a)
				
						normviol = (a*d - ((b + c)*(b + c)/4))/max(a*d,((b + c)*(b + c)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'a'));
					}
					else if(a < -interiortol && d < -interiortol && a*d - ((b + c)*(b + c)/4) > interiortol){ //this covers case (1b)
				
						normviol = (a*d - ((b + c)*(b + c)/4))/max(a*d,((b + c)*(b + c)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'b'));
					}
					else if( a > interiortol && d < -interiortol && -a*d - ((b - c)*(b - c)/4) > interiortol){ //this covers case (1c)
				
						normviol = (-a*d - ((b - c)*(b - c)/4))/max(-a*d,((b - c)*(b - c)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'c'));
					}
					else if(a < -interiortol && d > interiortol && -a*d - ((b - c)*(b - c)/4) > interiortol){ //this covers case (1d)
				
						normviol = (-a*d - ((b - c)*(b - c)/4))/max(-a*d,((b - c)*(b - c)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'd'));
					}

					if(b > interiortol && c > interiortol && b*c - ((a + d)*(a + d)/4) > interiortol){ //this covers case (1e)
				
						normviol = (b*c - ((a + d)*(a + d)/4))/max(b*c,((a + d)*(a + d)/4));

						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'e'));
					}
					else if(b < -interiortol && c < -interiortol && b*c - ((a + d)*(a + d)/4) > interiortol){ //this covers case (1f)
				
						normviol = (b*c - ((a + d)*(a + d)/4))/max(b*c,((a + d)*(a + d)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'f'));
					}
					else if(b > interiortol and c < -interiortol and -b*c - ((a - d)*(a - d)/4) > interiortol){ //this covers case (1g)
				
						normviol = (-b*c - ((a - d)*(a - d)/4))/max(-b*c,((a - d)*(a - d)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'g'));
					}
					else if(b < -interiortol and c > interiortol and -b*c - ((a - d)*(a - d)/4) > interiortol){ //this covers case (1h)
				
						normviol = (-b*c - ((a - d)*(a - d)/4))/max(-b*c,((a - d)*(a - d)/4));
						counter++;
						minor_violation_tuples.push_back(make_tuple(i,j,k,l, normviol, 'h'));
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
		
			MatrixXd D = avg_Dir;

			Vector2d linA(D(i,k),solX_matrix(i,k)); //a = lambda*D11 + X11
			Vector2d linB(D(i,l),solX_matrix(i,l)); //b = lambda*D12 + X12
			Vector2d linC(D(j,k),solX_matrix(j,k)); //c = lambda*D21 + X21
			Vector2d linD(D(j,l),solX_matrix(j,l)); //d = lambda*D22 + X22

			double steplength = std::numeric_limits<double>::infinity();

			if(set_type == 'a' || set_type == 'b'){

				

				if(abs(linA(0)) > eps && -linA(1)/linA(0) > eps)
					steplength = min(steplength, -linA(1)/linA(0)); //steplengths making a = 0
		
				if(abs(linD(0)) > eps && -linD(1)/linD(0) > eps)
					steplength = min(steplength, -linD(1)/linD(0)); //steplengths making d = 0

			

				Vector2d linSum = (linB + linC)/2;
				Vector3d quad = linearFormProd(linA,linD) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}
			else if(set_type == 'c' || set_type == 'd'){

				if(abs(linA(0)) > eps && -linA(1)/linA(0) > eps)
					steplength = min(steplength, -linA(1)/linA(0)); //steplengths making a = 0
		
				if(abs(linD(0)) > eps && -linD(1)/linD(0) > eps)
					steplength = min(steplength, -linD(1)/linD(0)); //steplengths making d = 0
			
				Vector2d linSum = (linB - linC)/2;
				Vector3d quad = - linearFormProd(linA,linD) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}
			else if(set_type == 'e' || set_type == 'f'){

				if(abs(linB(0)) > eps && -linB(1)/linB(0) > eps)
					steplength = min(steplength, -linB(1)/linB(0)); //steplengths making b = 0
		
				if(abs(linC(0)) > eps && -linC(1)/linC(0) > eps)
					steplength = min(steplength, -linC(1)/linC(0)); //steplengths making c = 0
			
				Vector2d linSum = (linA + linD)/2;
				Vector3d quad = linearFormProd(linB,linC) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));

			}
			else if(set_type == 'g' or set_type == 'h'){

				if(abs(linB(0)) > eps && -linB(1)/linB(0) > eps)
					steplength = min(steplength, -linB(1)/linB(0)); //steplengths making b = 0
		
				if(abs(linC(0)) > eps && -linC(1)/linC(0) > eps)
					steplength = min(steplength, -linC(1)/linC(0)); //steplengths making c = 0
			
				Vector2d linSum = (linA - linD)/2;
				Vector3d quad = - linearFormProd(linB,linC) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}

			get<4>(minor_violation_tuples[vind])= steplength;
		}
	}

	sort(begin(minor_violation_tuples), end(minor_violation_tuples),
			  [](tuple<int, int, int, int, double, char> const &t1, tuple<int, int, int, int, double, char> const &t2) 
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

		//find intersection points
		RowVectorXd Beta(N);
		Beta.setZero();
		
		for( int dir_num=0; dir_num < N; dir_num++){
			MatrixXd D = dirs_matrix[dir_num];

			Vector2d linA(D(i,k),solX_matrix(i,k)); //a = lambda*D11 + X11
			Vector2d linB(D(i,l),solX_matrix(i,l)); //b = lambda*D12 + X12
			Vector2d linC(D(j,k),solX_matrix(j,k)); //c = lambda*D21 + X21
			Vector2d linD(D(j,l),solX_matrix(j,l)); //d = lambda*D22 + X22

			double steplength = std::numeric_limits<double>::infinity();

			if(set_type == 'a' || set_type == 'b'){
				
				if(abs(linA(0)) > eps && -linA(1)/linA(0) > eps)
					steplength = min(steplength, -linA(1)/linA(0)); //steplengths making a = 0
		
				if(abs(linD(0)) > eps && -linD(1)/linD(0) > eps)
					steplength = min(steplength, -linD(1)/linD(0)); //steplengths making d = 0
			
				Vector2d linSum = (linB + linC)/2;
				Vector3d quad = linearFormProd(linA,linD) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}
			else if(set_type == 'c' || set_type == 'd'){

				if(abs(linA(0)) > eps && -linA(1)/linA(0) > eps)
					steplength = min(steplength, -linA(1)/linA(0)); //steplengths making a = 0
		
				if(abs(linD(0)) > eps && -linD(1)/linD(0) > eps)
					steplength = min(steplength, -linD(1)/linD(0)); //steplengths making d = 0
			
				Vector2d linSum = (linB - linC)/2;
				Vector3d quad = - linearFormProd(linA,linD) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}
			else if(set_type == 'e' || set_type == 'f'){

				if(abs(linB(0)) > eps && -linB(1)/linB(0) > eps)
					steplength = min(steplength, -linB(1)/linB(0)); //steplengths making b = 0
		
				if(abs(linC(0)) > eps && -linC(1)/linC(0) > eps)
					steplength = min(steplength, -linC(1)/linC(0)); //steplengths making c = 0
			
				Vector2d linSum = (linA + linD)/2;
				Vector3d quad = linearFormProd(linB,linC) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));

			}
			else if(set_type == 'g' or set_type == 'h'){

				if(abs(linB(0)) > eps && -linB(1)/linB(0) > eps)
					steplength = min(steplength, -linB(1)/linB(0)); //steplengths making b = 0
		
				if(abs(linC(0)) > eps && -linC(1)/linC(0) > eps)
					steplength = min(steplength, -linC(1)/linC(0)); //steplengths making c = 0
			
				Vector2d linSum = (linA - linD)/2;
				Vector3d quad = - linearFormProd(linB,linC) - linearFormProd(linSum,linSum);
				Vector2d roots = computeRoots(quad(0),quad(1),quad(2));

				if(roots(0) > eps)
					steplength = min(steplength, roots(0));

				if(roots(1) > eps)
					steplength = min(steplength, roots(1));
			}

			steplength *= (1-stepback);
			Beta(dir_num) = 1.0/steplength;    //make sure step keeps us in bounds
               
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


void minorcut(VectorXd solX, MatrixXd solX_matrix, vector<VectorXd> dirs, vector<MatrixXd> dirs_matrix, 
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic, int max_cuts,
		vector<RowVectorXd> *out_pi, vector<double> *out_pirhs, vector<double> *out_violation){

	//generate elementary minor cuts
	//to ensure separation, problem MUST have nonnegative diagonals

	double stepback=1E-8; //what fraction to step back to avoid infeas cuts
	double PDtol=1E-8;
	double interiortol = 1E-5; //what is considered to be in the interior of the 2x2PSD cone

	//find elementary violations
	int counter = 0;
	
	vector<tuple<int,int, double>> minor_violation_tuples;

	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n+1; j++){
			if(solX_matrix(i,i)*solX_matrix(j,j) - solX_matrix(i,j)*solX_matrix(i,j) > interiortol){
				
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
		MatrixXd lamDm;
		bool infflags[N];
		for(int i=0; i<N; i++)
			infflags[i] = false;

		for( int k=0; k < N; k++){
			MatrixXd D = dirs_matrix[k];

			Vector2d linA(D(ind1,ind1),solX_matrix(ind1,ind1)); //a = X11 + lambda*D11
			Vector2d linB(D(ind1,ind2),solX_matrix(ind1,ind2)); //b = X12 + lambda*D12
			Vector2d linC(D(ind2,ind2),solX_matrix(ind2,ind2)); //c = X22 + lambda*D22

			double steplength = std::numeric_limits<double>::infinity();

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
				lamDm = steplength*D;
			}
			else{
				infflags[k] = true;
			}

			steplength *= (1-stepback);
			
			Beta(k) = 1.0/steplength;
		}
               
		if(Beta.squaredNorm() < eps){
			cout << "Minor cut shows infeasible problem" << endl;
			exit(1);
		}

		//STRENGTHENING STARTS HERE
		
		if(strengthen){
			for(int i=0; i<N; i++){
				if(infflags[i]){
					MatrixXd D = dirs_matrix[i];

					/*
					tighten contained rays by rotation
					let C = constmat, D = dirmat
					want C+alpha*D to be exactly psd:
					(C11+alphaD11)(C22+alphaD22)=(C12+alphaD12)^2
					C11C22+alphaC11D22+alphaC22D11+alpha^2D11D22=C12^2+2alphaC12D12+alpha^2D12^2
					a=D11D22-D12^2
					b=C11D22+C22D11-2C12D12
					c=C11C22-C12^2
					ax^2+bx+c, want the nonpositive root closest to zero
					*/
					double a = D(ind1,ind1)*D(ind2,ind2) - D(ind1,ind2)*D(ind1,ind2);
					double b = - lamDm(ind1,ind1)*D(ind2,ind2) - lamDm(ind2,ind2)*D(ind1,ind1) + 2*lamDm(ind1,ind2)*D(ind1,ind2);
					double c = lamDm(ind1,ind1)*lamDm(ind2,ind2) - lamDm(ind1,ind2)*lamDm(ind1,ind2);

					Vector2d roots = computeRoots(a, b, c);

					if( min(roots(0), roots(1)) < 0 ) {
						double steplength = min(roots(0), roots(1));
						steplength *= (1+stepback); //note this has to be plus, you want this to be more negative

						Beta(i) = 1/steplength;
						//cout << "Direction strengthened" << endl;
						
		                    	}
					//else 
					//	cout << "Tilt error, positive steps, reverting " << min(roots(0), roots(1)) << endl ;//leave step at infinity

				}
			}
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
		int n, int N, int **Xtovec, MatrixXd Abasic, VectorXd bbasic,
		RowVectorXd *out_pi, double *out_pirhs, double *out_violation){
	
	double stepback = 1E-8;
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
            
			bool infflags[N];
			for(int i=0; i<N; i++)
				infflags[i] = false;

			bool testflag = true;
            		MatrixXd lamDm;

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
					/*
					if( roots == NULL ){
						cout << "Negative discriminant !! " << endl;
						*out_violation = -1;
						return;
					}
					*/
							
					if( min(roots(0), roots(1)) > 0 ) 
						steplength = min(roots(0), roots(1));  
					else
						steplength = max(roots(0), roots(1));
						
					steplength *= (1-stepback);
					
					MatrixXd testpt = solX_matrix + steplength*D;
                    
					double testr = myMatrixInnerProduct(testpt, C, n+1)*q/(Cfro*sqrt(Csqs));
					
					MatrixXd auxMatrix = testpt - myMatrixInnerProduct(testpt,C, n+1)*C/(Cfro*Cfro);
					double testd = auxMatrix.norm();
                    
					if(testr < testd){
						cout << "STEPLENGTH ERROR, REVERTING ---------" << endl;
						steplength = radius/D.norm();
					}	
					else if( (testr-testd)/testd > 0.01 ){
						cout << "big gap?" << endl;
						cout << testr << " vs " << testd <<endl;
						if(radius/D.norm() > steplength){
							cout << "steplength is smaller than ballstep, steplength: " << steplength << endl;
							cout << "ballstep: " << radius/D.norm();
							steplength = radius/D.norm();
						}
					}
					
					lamDm = steplength*D; //for strengthening
					testflag=0;
				}	
				else{
					//INTERSECTION AT INFINITY
					steplength = std::numeric_limits<double>::infinity();
					infflags[i] = true;
				}
					
				Beta(i) = 1/steplength;
			}

			//strengthen infinite steps
			if(strengthen){
				double innerLC = myMatrixInnerProduct(lamDm,C, n+1);
				MatrixXd Z1 = Cfro*lamDm - innerLC*C/Cfro;
		    
				for(int i=0; i<N; i++){
					if(infflags[i]){
						MatrixXd D = dirs_matrix[i];

						double innerDC = myMatrixInnerProduct(D,C,n+1);
						MatrixXd Z2 = innerDC*C/Cfro-Cfro*D;
						double quada = (q*q*innerDC*innerDC)/Csqs - Z2.squaredNorm();
						double quadb = -2*(q*q)*innerLC*innerDC/Csqs - 2*myMatrixInnerProduct(Z1,Z2,n+1);
						double quadc = (q*q*innerLC*innerLC)/Csqs - Z1.squaredNorm();

						double steplength = 0;                    

						Vector2d roots = computeRoots(quada, quadb, quadc);
						/*
						if( roots == NULL){
							cout << "Negative discriminant !! " << endl;
							*out_violation = -1;
							return;
						}
						*/
						if( min(roots(0), roots(1)) < 0 ) {
							steplength = min(roots(0), roots(1));
							steplength *= (1+stepback); //note this has to be plus, you want this to be more negative
							MatrixXd testdir = lamDm - steplength*D;
						
							double rtest = myMatrixInnerProduct(testdir,C,n+1)*q/(Cfro*sqrt(Cfro*Cfro - q*q));
						
							MatrixXd auxMatrix = testdir - myMatrixInnerProduct(testdir,C, n+1)*C/(Cfro*Cfro);
							double dtest = auxMatrix.norm();

							if(dtest <= rtest)
								Beta(i) = 1/steplength;
							else
								cout << "Tilt direction not in recession cone, reverting" << endl; //leave step at infinity
		                    		}
						else 
							cout << "Tilt error, positive steps, reverting" << endl ;//leave step at infinity
					}
				}
			}
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

	double stepback = 1E-8;

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
