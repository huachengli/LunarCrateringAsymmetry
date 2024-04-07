//
// Created by huachengli on 2020/8/15.
//

#ifndef GAUSSIANINT3_INTEGRATE_BASE_H
#define GAUSSIANINT3_INTEGRATE_BASE_H
#include <math.h>
#include <omp.h>
#include <stdio.h>
#define getbit(x,y)   ((x) >> (y)&1)
double BaseIntGi5(double * Xa,double * Xb,int dim,double (*f)(double *,int,double *,int),double * uav,int uac);
double IntGi5(double ** X,int *nX,int dim,double (*f)(double *,int,double *,int),double * uav,int uac);
#endif //GAUSSIANINT3_INTEGRATE_BASE_H
