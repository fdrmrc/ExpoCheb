//============================================================================
// Name        : ExpoChebychev.cpp
//============================================================================

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <Eigen>
#include "mpi.h"
//#include "petsc.h"

/* my headers */
#include "ChebCoeffs.h"
#include "integrator.h"

using namespace std;
using namespace Eigen;

MatrixXd build_matrix(int n,double h);
int main(int argc, char* argv[]) {
	cout << "\t Chebychev rational approximation method \t \n" << endl;

	complex<double> alpha[p];
	alpha[0] = 1.8321743782540412751 * 1e-14;
	alpha[1] = -7.1542880635890672853 * 1e-5 + i*1.4361043349541300111 * 1e-4;
	alpha[2] =  +9.4390253107361688779* 1e-3  - i*1.7184791958483017511 * 1e-2;
	alpha[3] = -3.7636003878226968717 *1e-1 + i*3.3518347029450104214 * 1e-1;
	alpha[4] = -2.3498232091082701191* 1e1 - i*5.8083591297142074004 * 1e0;
	alpha[5] = +4.6933274488831293047 * 1e1 + i*4.5643649768827760791 * 1e1;
	alpha[6] = -2.7875161940145646468 * 1e1 - i*1.0214733999056451434 * 1e2;
	alpha[7] = +4.8071120988325088907 * 1e0 - i*1.3209793837428723881 * 1e0;


	complex<double> theta[p-1];
	theta[0] = - 8.8977731864688888199 * 1e0 + i*1.6630982619902085304 * 1e1;
	theta[1] = - 3.7032750494234480603 * 1e0 + i*1.3656371871483268171 * 1e1;
	theta[2] = - 0.2087586382501301251 * 1e0 + i*1.0991260561901260913 * 1e1;
	theta[3] = + 3.9933697105785685194 * 1e0 + i*6.0048316422350373178 * 1e0;
	theta[4] = + 5.0893450605806245066 * 1e0 + i*3.5888240290270065102 * 1e0;
	theta[5] = + 5.6231425727459771248 * 1e0 + i*1.1940690463439669766 * 1e0;
	theta[6] = + 2.2697838292311127097 * 1e0 + i*8.4617379730402214019 * 1e0;

	//test
	int dim = 10;
	double h = 0.1;
	MatrixXd B = build_matrix(dim,h); //laplacian matrix
	expoCheb eAv(B, dim,h);
	VectorXd v = eAv.build_vec();
	//eAv.stamp_mat(B); //visualize matrix


	VectorXd eLv = eAv.actexp(v,alpha,theta);
	eAv.stamp_vec(eLv);
        
        
        integrator myintegr(1.0,100,B,2.0*v,v);
        VectorXd a = myintegr.stepper(v,alpha,theta);
        cout << a << endl;
        

	return 0;
}


MatrixXd build_matrix(int n,double h){
	//build the laplacian matrix
	MatrixXd M(n,n);
	double h2 = h*h;
	for (int i=0;i<n;i++){
		M(i,i) = -2.0/h2; //diagonal
	}

	for (int j=1;j<n-1;j++){
		M(j,j-1) = 1.0/h2;
		M(j,j+1) = 1.0/h2;
	}
	M(0,1) = 1.0/h2;
	M(n-1,n-2) = 1.0/h2;
	return M;
}