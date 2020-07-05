/*
 * ChebCoeffs.h
 *
 *      Author: marco
 */

#ifndef CHEBCOEFFS_H_
#define CHEBCOEFFS_H_

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <Eigen>
#include "mpi.h"
//#include "ChebCoeffs.h"

using namespace std;
using namespace Eigen;


const complex<double> i(0.0,1.0);

extern const int p = 8;
extern complex<double> alpha[];
extern complex<double> theta[];


class expoCheb{
private:
	MatrixXd A;
	int n;
	double h;

public:
	expoCheb(MatrixXd matrix,int dim,double step){
			A = matrix;
			n = dim;
			h = step;
}


void refine_step(double ref, double &h){
	h/=ref;
}


VectorXd build_vec(){
	//VectorXcd vec = VectorXcd::Random(n);
	VectorXd vec = VectorXd::Ones(n);
	return vec;
}


void stamp_mat(MatrixXd M){
	cout << "\n \t MATRIX is \t \n"<< M << endl;
}

void stamp_vec(VectorXd vec){
	cout << "\n \t SOLUTION is \t \n"<< vec << endl;
}

VectorXcd lin_solver(MatrixXcd M,VectorXcd v){
	//solve the linear system Ax=v with Eigen Dense routines
	VectorXcd sol = M.fullPivLu().solve(v);
	return sol;
}

VectorXd actexp(VectorXcd v,complex<double> alpha[],complex<double> theta[]){
	VectorXcd sum;
	sum = alpha[0]*v;
	MatrixXd I = MatrixXd::Identity(n,n);
	for (int j=0;j<p-1;j++){
		//cout << j << endl;
		sum += 2.0*alpha[j+1]*lin_solver(A-theta[j]*I,v);
	}

	return sum.real();
}



};

#endif /* CHEBCOEFFS_H_ */