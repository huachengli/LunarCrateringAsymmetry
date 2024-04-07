//
// Created by lhc on 2/4/2024.
//
// collect some model function from the original main files

#ifndef GAUSSIANINT3_INTEGRATE_MODEL_H
#define GAUSSIANINT3_INTEGRATE_MODEL_H

#include<math.h>
double Ap(double plon,double plat);
double f1(double hx,double ga);
double g1(double hx,double ga);
double f2(double hx,double ga);
double g2(double hx,double ga);
double flux(double *X,int dim,double * uav,int uac);
double calvp(double * res,double *X,double * uav,int uac);
double calnm(double * res,double *X,double * uav,int uac);
double calvm(double * res,double *X,double * uav,int uac);
double nspeed(double *X,int dim,double * uav,int uac);
#endif //GAUSSIANINT3_INTEGRATE_MODEL_H
