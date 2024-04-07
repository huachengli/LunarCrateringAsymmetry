//
// Created by lhc on 2/4/2024.
//

#include "integrate_model.h"

double Ap(double plon,double plat)
{
    // the inclination of asteroids (relative to the ecliptic) that encounter with the Moon
    static double p[100] = {
            0.10789000, 0.10789193, 0.10728255, 0.10652964, 0.10541574, 0.10443136,  \
            0.10341175, 0.10239125, 0.10150454, 0.10081894, 0.10029208, 0.09988673,  \
            0.09958001, 0.09931634, 0.09890344, 0.09796903, 0.09656588, 0.09532840,  \
            0.09405840, 0.09278840, 0.09113516, 0.08885395, 0.08676451, 0.08475458,  \
            0.08273222, 0.08102610, 0.07987958, 0.07866554, 0.07754241, 0.07674952,  \
            0.07544289, 0.07369997, 0.07253604, 0.07097247, 0.06929820, 0.06769347,  \
            0.06608873, 0.06464497, 0.06295300, 0.06144448, 0.05994834, 0.05867226,  \
            0.05803839, 0.05749012, 0.05694715, 0.05646707, 0.05604520, 0.05582724,  \
            0.05581317, 0.05573929, 0.05560017, 0.05559969, 0.05527568, 0.05423406,  \
            0.05274399, 0.05147409, 0.04995946, 0.04843502, 0.04708107, 0.04602317,  \
            0.04506717, 0.04436267, 0.04312870, 0.04190157, 0.04098117, 0.03981434,  \
            0.03825078, 0.03719285, 0.03613495, 0.03512069, 0.03421890, 0.03337297,  \
            0.03252707, 0.03189328, 0.03099840, 0.02945118, 0.02829548, 0.02705810,  \
            0.02553128, 0.02407784, 0.02280783, 0.02163223, 0.02063493, 0.01984211,  \
            0.01884372, 0.01780374, 0.01666099, 0.01539429, 0.01398298, 0.01269945,  \
            0.01142951, 0.01015957, 0.00877714, 0.00739953, 0.00619598, 0.00511243,  \
            0.00402999, 0.00284015, 0.00261107, 0.00000000
    };
    static double dlat = M_PI_2/99.0;
    double lat = fabs(plat);
    int k = (int) floor(lat/dlat);
    if(k>=0 && k<= 98)
    {
        return (p[k] + (lat/dlat - k)*(p[k+1] - p[k]))/M_PI * 0.5 /0.174787 ;
        // 0.174787 is the result of trapz p[100]
        // to make int[p(x,y),{x,y}] = 1
    } else
        return 0;
}

double test1(double * X,int dim,double * uav,int uac)
{
    double plon = X[0];
    double plat = X[1];
    return Ap(plon,plat);
//    return 0.5*cos(plat)/M_PI;
}

double f1(double hx,double ga)
{
    if(hx < 0)
        return 0;
    else
        return hx;
}

double f2(double hx,double ga)
{
    if(hx < -ga/(2 + ga))
        return 0;
    else
        return hx + (0.5 - 2*hx)*ga;
}

double g1(double hx,double ga)
{
    return hx;
}

double g2(double hx,double ga)
{
    if(hx > 0)
        return hx + ga/2.0;
    else
        return -hx - ga/2.0;
}

double calvp(double * res,double *X,double * uav,int uac)
{
    // calculate the unit vector
    // Rz(plon - o2)*Rx(pi/2 - plat)*ez

    double plon = X[0];
    double plat = X[1];
    double o2   = X[2] + M_PI;
    // o2 = o1 + M_PI
    // Vp_x;Vp_y;Vp_z;
    res[0] =   cos(plat) * sin(plon - o2);
    res[1] = - cos(plat) * cos(plon - o2);
    res[2] =   sin(plat);
    return res[2];
}

double calnm(double * res,double *X,double * uav,int uac)
{
    // calculate the unit vector
    // Rx(i2)*Rz(pi/2+M+M0+lon)*Rx(pi/2 - lat)*ez

    double o1   = X[2];
    double o2   = X[2] + M_PI;
    double o3   = X[3];
    double E    = X[4];
    double e    = uav[0];
    double i1   = uav[1];
    double i2   = uav[2];
    double M    = E - e * sin(E);
    double lon  = uav[4];
    double lat  = uav[5];
    // considering the different type orbit resonance
    // M = 1.5*M;

    double tmpA = cos(i1)*sin(o1-o2)*sin(o3) - cos(o1 - o2)*cos(o3);
//    double tmpB = -cos(i1)*cos(i2)*cos(o1-o2)*sin(o3) - cos(i2)*sin(o1-o2)*cos(o3) + sin(i1)*sin(i2)*sin(o3);
//    correct for rotation sin(t)
    double tmpB = -cos(i1)*cos(i2)*cos(o1-o2)*sin(o3) - cos(i2)*sin(o1-o2)*cos(o3) - sin(i1)*sin(i2)*sin(o3);
    double cM0  = tmpA/sqrt(tmpA*tmpA + tmpB*tmpB);
    double sM0  = tmpB/sqrt(tmpA*tmpA + tmpB*tmpB);
    double cMM0 = cos(M + lon)*cM0 - sin(M + lon)*sM0;
    double sMM0 = sin(M + lon)*cM0 + cos(M + lon)*sM0;
    // Nm_x;Nm_y;Nm_z
    res[0] =  cMM0 * cos(lat);
    res[1] =  cos(i2)*sMM0*cos(lat) - sin(i2)*sin(lat);
    //  res[2] = -sin(i2)*sMM0*cos(lat) + cos(i2)*sin(lat);
    //  correction for rotation sin(t)
    res[2] = sin(i2)*sMM0*cos(lat) + cos(i2)*sin(lat);
    return res[2];
}

double calvm(double * res,double *X,double * uav,int uac)
{
    // calculate a unit vector
    // Rz(o1-o2)*Rx(i1)*Rz(o3)*v
    // v = [ vm*sin(f),vm*(e + cos(f)),0]^T

    double o1   = X[2];
    double o2   = X[2] + M_PI;
    double o3   = X[3];
    double E    = X[4];
    double e    = uav[0];
    double i1   = uav[1];
    double cfm   = (cos(E) - e)/(1 - e*cos(E));
    double sfm   = sqrt(1 - e*e) * sin(E) / (1 - e*cos(E));
    double eta  = uav[3];

    // z(o1-o2)*Rx(i1)*Rz(o3)*ex
    double vxx = cos(o1-o2)*cos(o3) - cos(i1)*sin(o1-o2)*sin(o3);
    double vxy = sin(o1-o2)*cos(o3) + cos(i1)*cos(o1-o2)*sin(o3);
//    double vxz = -sin(i1)*sin(o3);
//    correction for rotation sin(t)
    double vxz =  sin(i1)*sin(o3);

    double vyx = -cos(i1)*sin(o1-o2)*cos(o3) - cos(o1-o2)*sin(o3);
    double vyy =  cos(i1)*cos(o1-o2)*cos(o3) - sin(o1-o2)*sin(o3);
//    double vyz = -sin(i1)*cos(o3);
//    correction for rotation sin(t)
    double vyz = sin(i1)*cos(o3);

    double vmx = - eta/sqrt(1 - e*e) * sfm;
    double vmy = eta/sqrt(1 - e*e) * (e + cfm);
    res[0] = vmx*vxx + vmy*vyx;
    res[1] = vmx*vxy + vmy*vyy;
    res[2] = vmx*vxz + vmy*vyz;
    return res[2];
}

double flux(double *X,int dim,double * uav,int uac)
{
    double plon = X[0];
    double plat = X[1];
    double E    = X[4];
    double e    = uav[0];
    double gamma0 = uav[6];
    double vp[3];
    double vm[3];
    double nm[3];
    calnm(nm,X,uav,uac);
    calvm(vm,X,uav,uac);
    calvp(vp,X,uav,uac);
    double Nvpvm = (vp[0] - vm[0])*(vp[0]-vm[0]) + (vp[1] - vm[1])*(vp[1]-vm[1]) + (vp[2] - vm[2])*(vp[2]-vm[2]);
    double hatx = (vm[0]-vp[0])*nm[0] + (vm[1]-vp[1])*nm[1] + (vm[2]-vp[2])*nm[2];
    hatx = hatx/sqrt(Nvpvm);
    double gamma = gamma0/Nvpvm;
    // flux
    // C * |vp-vm| * p(plon,plat) * f(hatx) * dM
    // dM = (1 - e*cos(E))dE
    return sqrt(Nvpvm)*f2(hatx,gamma)*Ap(plon,plat)*(1 - e*cos(E));
}

double nspeed(double *X,int dim,double * uav,int uac)
{
    double plon = X[0];
    double plat = X[1];
    double E    = X[4];
    double e    = uav[0];
    double gamma0 = uav[6];
    double vp[3];
    double vm[3];
    double nm[3];
    calnm(nm,X,uav,uac);
    calvm(vm,X,uav,uac);
    calvp(vp,X,uav,uac);
    double Nvpvm = (vp[0] - vm[0])*(vp[0]-vm[0]) + (vp[1] - vm[1])*(vp[1]-vm[1]) + (vp[2] - vm[2])*(vp[2]-vm[2]);
    double hatx = (vm[0]-vp[0])*nm[0] + (vm[1]-vp[1])*nm[1] + (vm[2]-vp[2])*nm[2];
    hatx = hatx/sqrt(Nvpvm);
    double gamma = gamma0/Nvpvm;
    // flux
    // C * |vp-vm| * p(plon,plat) * f(hatx) * dM
    // dM = (1 - e*cos(E))dE
    return Nvpvm*f2(hatx,gamma)*Ap(plon,plat)*(1 - e*cos(E)) * g2(hatx,gamma);
}
