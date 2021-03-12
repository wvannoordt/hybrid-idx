#ifndef MMS_H
#define MMS_H

#define PI (3*3.14159265359)

//pi/4
#define XYROT 0.7

#define xrot(x, y) ((x)*cos(XYROT) - (y)*sin(XYROT))
#define yrot(x, y) ((x)*sin(XYROT) + (y)*cos(XYROT))

#define pres_mms(ar, x, y, z) {ar[0] = 4.0 + cos(PI*(xrot(x,y)))*cos(PI*(0.0));\
ar[1] = -PI*sin(PI*(xrot(x, y)))*cos(XYROT)*cos(PI*(0.0));\
ar[2] =  PI*sin(PI*(xrot(x, y)))*sin(XYROT)*cos(PI*(0.0));\
ar[3] = 0.0;\
}

#define dens_mms(ar, x, y, z) {ar[0] = 4.0 + cos(PI*(xrot(x,y)))*cos(PI*(0.0));\
ar[1] = -PI*sin(PI*(xrot(x, y)))*cos(XYROT)*cos(PI*(0.0));\
ar[2] =  PI*sin(PI*(xrot(x, y)))*sin(XYROT)*cos(PI*(0.0));\
ar[3] = 0.0;\
}

#define uvel_mms(ar, x, y, z) {ar[0] = 4.0 + cos(PI*(xrot(x,y)))*cos(PI*(0.0));\
ar[1] = -PI*sin(PI*(xrot(x, y)))*cos(XYROT)*cos(PI*(0.0));\
ar[2] =  PI*sin(PI*(xrot(x, y)))*sin(XYROT)*cos(PI*(0.0));\
ar[3] = 0.0;\
}

#define vvel_mms(ar, x, y, z) {ar[0] = 4.0 + cos(PI*(xrot(x,y)))*cos(PI*(0.0));\
ar[1] = -PI*sin(PI*(xrot(x, y)))*cos(XYROT)*cos(PI*(0.0));\
ar[2] =  PI*sin(PI*(xrot(x, y)))*sin(XYROT)*cos(PI*(0.0));\
ar[3] = 0.0;\
}

#define wvel_mms(ar, x, y, z) {ar[0] = 4.0 + cos(PI*(xrot(x,y)))*cos(PI*(0.0));\
ar[1] = -PI*sin(PI*(xrot(x, y)))*cos(XYROT)*cos(PI*(0.0));\
ar[2] =  PI*sin(PI*(xrot(x, y)))*sin(XYROT)*cos(PI*(0.0));\
ar[3] = 0.0;\
}

#define cont_rhs_mms(p, rho, u, v, w) (u[0]*rho[1]+v[0]*rho[2]+w[0]*rho[3] + rho[0]*(u[1] + v[2] + w[3]))
    
#define engy_rhs_mms(p, rho, u, v, w, e) (-u[0]*p[1]-p[0]*u[1]-v[0]*p[2]-p[0]*v[2]-w[0]*p[3]-p[0]*w[3]\
    - u[0]*e[0]*rho[1] - rho[0]*(e[0]*u[1] + u[0]*e[1]) - e[0]*v[0]*rho[2] - rho[0]*(e[0]*v[2] + v[0]*e[2]) - e[0]*w[0]*rho[3] - rho[0]*(e[0]*w[3] + w[0]*e[3]))

#define momx_rhs_mms(p, rho, u, v, w) (-p[1] - u[0]*u[0]*rho[1] - rho[0]*(u[0]*u[1] + u[0]*u[1]) - u[0]*v[0]*rho[2] - rho[0]*(u[0]*v[2] + v[0]*u[2]) - u[0]*w[0]*rho[3] - rho[0]*(u[0]*w[3] + w[0]*u[3]))
#define momy_rhs_mms(p, rho, u, v, w) (-p[2] - v[0]*u[0]*rho[1] - rho[0]*(v[0]*u[1] + u[0]*v[1]) - v[0]*v[0]*rho[2] - rho[0]*(v[0]*v[2] + v[0]*v[2]) - v[0]*w[0]*rho[3] - rho[0]*(v[0]*w[3] + w[0]*v[3]))
#define momz_rhs_mms(p, rho, u, v, w) (-p[3] - w[0]*u[0]*rho[1] - rho[0]*(w[0]*u[1] + u[0]*w[1]) - w[0]*v[0]*rho[2] - rho[0]*(w[0]*v[2] + v[0]*w[2]) - w[0]*w[0]*rho[3] - rho[0]*(w[0]*w[3] + w[0]*w[3]))

#define sqr(x) (x)*(x)

#define d_abs(x) ((((x)<0)?(-(x)):((x))))

#endif