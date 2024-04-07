#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "integrate_base.h"
#include "integrate_model.h"

int main()
{
    struct timeval start,end;
    gettimeofday(&start, NULL );

    // build integrate grid for Int5
    double c0 = 8*pow(M_PI,3);
    int nX[5] = {31,31,11,11,11};
    double Xa[5] = {-M_PI,-M_PI_2,0,0,0};
    double Xb[5] = {M_PI,M_PI_2,2*M_PI,2*M_PI,2*M_PI};
    double ** X = (double **)malloc(sizeof(double*)*5);

    for(int k=0;k<5;k++)
    {
        X[k] = (double *) malloc(sizeof(double)*nX[k]);
    }
    for(int k=0;k<5;k++)
    {
        for(int i=0;i<nX[k];i++)
        {
            X[k][i] = Xa[k] + (Xb[k] - Xa[k])/(nX[k]-1.0) * i;
        }
    }
    double deg2rad = M_PI/180.0;

    // build uniform grid on (lon,lat)
    double x[2] = {-M_PI,M_PI};
    double y[2] = {-M_PI_2,M_PI_2};
    int nx = 90, ny = 45;
    double * xs = (double *)malloc(sizeof(double)*nx);
    double * yt = (double *)malloc(sizeof(double)*ny);
    for(int i=0;i<nx;i++) xs[i] = x[0] + (x[1]-x[0]) * i /(nx - 1);
    for(int j=0;j<ny;j++) yt[j] = y[0] + (y[1]-y[0]) * j /(ny -1 );

    // some physical parameters for Earth-Moon system
    double Re = 6.371e6;
    double Rm = 1.73740e6;
    double GMe = 5.972e24 * 6.674e-11;
    double GMm = 7.342e22 * 6.674e-11;
    double vm0 = sqrt(GMe/Re);
    double vp = 19.0*1000.0;
    double vm = 1.022*1000.0;
    double Gamma0 = 2*GMm/Rm / (vp*vp);
    double uav[7] = {0.0549,5.145*deg2rad,1.535*deg2rad,vm/vp,-M_PI_2,0,Gamma0};
    int uac = 6;
    // some parameters that need attached to nspeed/flux
    // e ; i1 ; i2 ; eta ; lon ; lat ; gamma0
    // eta = vm at perige / vp
    // gamma = 2GMm/RU^2

    double ** S;
    int Ssize = nx*ny;
    S = (double **) malloc(sizeof(double*)*Ssize);
    for(int i=0;i<Ssize;i++)
    {
        S[i] = (double *) malloc(sizeof(double)*4);
    }

    for(int ix=0;ix<nx;ix++)
    {
        for(int iy=0;iy<ny;iy++)
        {
            int Index = ix*ny + iy;
            S[Index][0] = xs[ix];
            S[Index][1] = yt[iy];
        }
    }

    // calculate speed/flux distribution
    fprintf(stdout,"calculating speed/flux distribution on (%d,%d) grid ...\n",nx,ny);
    #pragma omp parallel for num_threads(16)
    for(int i=0;i<Ssize;i++)
    {
        double louav[7] = {0.0549,5.145*deg2rad,1.535*deg2rad,vm/vp,S[i][0],S[i][1],Gamma0};
        S[i][2] = IntGi5(X,nX,5,flux,louav,uac)/c0;
        S[i][3] = IntGi5(X,nX,5,nspeed,louav,uac)/c0;
    }

    // write output;
    char fname[40];
    sprintf(fname,"lon_lat_am60.csv");
    fprintf(stdout,"writing %s ....\n",fname);
    FILE * fp = fopen(fname,"w");
    assert(fp!=NULL);

    for(int i=0;i<Ssize;i++)
    {
        for(int k=0;k<4;k++)
            fprintf(fp,"%4.5f, ",S[i][k]);
        fprintf(fp,"\n");
    }
    fclose(fp);

    for(int k=0;k<5;k++) free(X[k]);
    free(X);
    for(int k=0;k<Ssize;k++) free(S[k]);
    free(S);
    free(xs);free(yt);

    gettimeofday(&end, NULL );
    double ds = (end.tv_sec-start.tv_sec) + (end.tv_usec-start.tv_usec)/1000000.0;
    fprintf(stdout,"end Total time: %6.5f sec",ds);
    return 0;
}


