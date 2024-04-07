//
// Created by huachengli on 2020/8/15.
//
// edited 2024/04/02

#include "integrate_base.h"

double BaseIntGi5(double * Xa,double * Xb,int dim,double (*f)(double *,int,double *,int),double * uav,int uac)
{
    double gpt= 1.0/sqrt(3);
    if (5!=dim)
    {
        fprintf(stdout,"the dimension is %d in BaseIntGi5",dim);
        return 0;
    }

    double midX[5] = {0};
    double dX[5]   = {0};
    double res     = 0.0;
    double dS      = 1.0;
    for(int i=0;i<5;i++)
    {
        midX[i] = 0.5*(Xb[i] + Xa[i]);
        dX[i]   = 0.5*(Xb[i] - Xa[i]);
        dS      = dS*dX[i];
    }
    double tmpX[5] = {0};
    for(int k=0;k<32;k++)
    {
        for(int i=0;i<5;i++)
        {
            tmpX[i] = (2.0*getbit(k,i) -1) * gpt;
            tmpX[i] = midX[i] + dX[i] * tmpX[i];
        }
        res = res + f(tmpX,5,uav,uac);
    }
    return res*dS;
}

double IntGi5(double ** X,int *nX,int dim,double (*f)(double *,int,double *,int),double * uav,int uac)
{
    double res =0;
    int i[5] = {0};

    if (5!=dim)
    {
        fprintf(stdout,"the dimension is %d in BaseIntGi5",dim);
        return 0;
    }

    for(int k=0;k<nX[0]*nX[1];k++)
    {
        int i0 = k/nX[1];
        int i1 = k%nX[1];
        double tmpa[5];
        double tmpb[5];
        tmpa[0] = X[0][i0];
        tmpb[0] = X[0][i0+1];
        tmpa[1] = X[1][i1];
        tmpb[1] = X[1][i1+1];

        if(nX[0]-1 != i0 && nX[1]-1 != i1)
        {
            for(int i2=0;i2<nX[2]-1;i2++)
            {
                tmpa[2] = X[2][i2];
                tmpb[2] = X[2][i2+1];
                for(int i3=0;i3<nX[3]-1;i3++)
                {
                    tmpa[3] = X[3][i3];
                    tmpb[3] = X[3][i3+1];
                    for(int i4=0;i4<nX[4]-1;i4++)
                    {
                        tmpa[4] = X[4][i4];
                        tmpb[4] = X[4][i4+1];
                        res += BaseIntGi5(tmpa,tmpb,dim,f,uav,uac);
                    }
                }
            }
        }
    }
    return res;
}
