#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <cmath>
#include "ChebCoeffs.h"
#include <Eigen>


class integrator{
private:
    double T;
    int ts; 
    MatrixXd A;
    VectorXd g;
    VectorXd y0;
 
public:
    integrator(double fin_time, int time_steps, MatrixXd rhs, VectorXd frhs, VectorXd in_data){
        T = fin_time;
        ts = time_steps;
        A = rhs;
        g = frhs;
        y0 = in_data;
    }

double get_step(){
    double k = T/ts;
    return k;
}
    
VectorXd stepper(VectorXd Un,complex<double> alpha[],complex<double> theta[]){
    VectorXd ret,Aun;
    double k = get_step();
    Aun = A*Un;
    expoCheb eAv(0.5*k*A, Un.rows(), 0.1);
    VectorXd xi1 = Un;
    VectorXd xi2 = Un  + 0.5*k*eAv.actexp(g + A*Un,alpha,theta);
    ret = Un + k*eAv.actexp(g + Aun,alpha,theta);
    return ret;
}



};

#endif /* INTEGRATOR_H */