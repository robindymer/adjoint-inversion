
#include "mex.h"
#include <omp.h>

/* The computational routines */
void D2_2nd(double *u_x, double *u, double *b, double dx, int Nx, int Ny)
{
    int i, j, bp, idx, str;
    double dxm = 1.0/(dx*dx);

    double wm1, w0, wp1;
    bp = 2;

    str = Ny;
    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            idx = j*Ny + i;

            wm1 = 0.5*(b[idx - str] + b[idx]);
            w0 = -0.5*(b[idx - str] + b[idx + str]) - b[idx];
            wp1 = 0.5*(b[idx + str] + b[idx]);

            u_x[idx] = dxm*(wm1*u[idx - str] + w0*u[idx] + wp1*u[idx + str]);
        }
    }
}

void D2_4th(double *u_x, double *u, double *b, double dx, int Nx, int Ny)
{
    int i, j, bp, str, str2;
    int im2, im1, i0, ip1, ip2;
    double dxm = 1.0/(dx*dx);

    double wm2, wm1, w0, wp1, wp2;
    bp = 6;

    double c1 = 1.0/8.0;
    double c2 = 1.0/6.0;
    double c3 = 1.0/24.0;
    double c4 = 1.0/1.2;

    str = Ny;
    str2 = 2*str;
    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            i0 = j*Ny + i;
            im2 = i0-str2;
            im1 = i0-str;
            ip1 = i0+str;
            ip2 = i0+str2;

            wm2 = -(b[im2] + b[i0])*c1 + b[im1]*c2;
            wm1 = (b[im2] + b[ip1])*c2 + (b[im1] + b[i0])*0.5;
            w0 = -(b[im2] + b[ip2])*c3 -(b[im1] + b[ip1])*c4 - 0.75*b[i0];
            wp1 = (b[im1] + b[ip2])*c2 + (b[i0] + b[ip1])*0.5;
            wp2 = -(b[i0] + b[ip2])*c1 + b[ip1]*c2;

            u_x[i0] = dxm*(wm2*u[im2] + wm1*u[im1] + w0*u[i0] + wp1*u[ip1] + wp2*u[ip2]);
        }
    }
}

void D2_6th(double *out, double *u, double *bx, double *by, double dx, double dy, int Nx, int Ny)
{
    int i, j, bp, str, str2, str3;
    int im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/(dx*dx);
    double dym = 1.0/(dy*dy);

    double wm3, wm2, wm1, w0, wp1, wp2, wp3;
    bp = 9;

    double c1 = -dxm*1.0/40.0;
    double c2 = -dxm*0.11e2/0.360e3;
    double c3 = -dxm*1.0/0.20e2;
    double c4 = -dxm*0.3e1/0.10e2;
    double c5 = 0.7e1*c1;
    double c6 = 0.17e2*c1;
    double c7 = -dxm*1.0/0.180e3;
    double c8 = -dxm*0.101e3/0.180e3;
    double c9 = -dxm*1.0/0.8e1;
    double c10 = 0.19e2*c3;

    double d1 = -dym*1.0/40.0;
    double d2 = -dym*0.11e2/0.360e3;
    double d3 = -dym*1.0/0.20e2;
    double d4 = -dym*0.3e1/0.10e2;
    double d5 = 0.7e1*d1;
    double d6 = 0.17e2*d1;
    double d7 = -dym*1.0/0.180e3;
    double d8 = -dym*0.101e3/0.180e3;
    double d9 = -dym*1.0/0.8e1;
    double d10 = 0.19e2*d3;


    #pragma omp parallel for private(i,j,im3,im2,im1,i0,ip1,ip2,ip3,wm3,wm2,wm1,w0,wp1,wp2,wp3)
    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            i0 = j*Ny + i;

            im3 = i0-3*Ny;
            im2 = i0-2*Ny;
            im1 = i0-Ny;
            ip1 = i0+Ny;
            ip2 = i0+2*Ny;
            ip3 = i0+3*Ny;

            wm3 = -c2*(bx[im3]+bx[i0]) + c1*(bx[im2]+bx[im1]);
            wm2 = c3*(bx[im3]+bx[ip1]) + c5*(bx[im2]+bx[i0]) - c4*bx[im1];
            wm1 = -c1*(bx[im3]+bx[ip2]) - c4*(bx[im2]+bx[ip1]) - c6*(bx[im1]+bx[i0]);
            w0  =  c7*(bx[im3]+bx[ip3]) + c9*(bx[im2]+bx[ip2]) + c10*(bx[im1]+bx[ip1]) + c8*bx[i0];
            wp1 = -c1*(bx[im2]+bx[ip3]) - c4*(bx[im1]+bx[ip2]) - c6*(bx[i0]+bx[ip1]);
            wp2 =  c3*(bx[im1]+bx[ip3]) + c5*(bx[i0]+bx[ip2]) - c4*bx[ip1];
            wp3 =  -c2*(bx[i0]+bx[ip3]) + c1*(bx[ip1]+bx[ip2]);

            out[i0] = wm3*u[im3] + wm2*u[im2] + wm1*u[im1] + w0*u[i0] + wp1*u[ip1] + wp2*u[ip2] + wp3*u[ip3];

            im3 = i0-3;
            im2 = i0-2;
            im1 = i0-1;
            ip1 = i0+1;
            ip2 = i0+2;
            ip3 = i0+3;

            wm3 = -d2*(by[im3]+by[i0]) + d1*(by[im2]+by[im1]);
            wm2 = d3*(by[im3]+by[ip1]) + d5*(by[im2]+by[i0]) - d4*by[im1];
            wm1 = -d1*(by[im3]+by[ip2]) - d4*(by[im2]+by[ip1]) - d6*(by[im1]+by[i0]);
            w0  =  d7*(by[im3]+by[ip3]) + d9*(by[im2]+by[ip2]) + d10*(by[im1]+by[ip1]) + d8*by[i0];
            wp1 = -d1*(by[im2]+by[ip3]) - d4*(by[im1]+by[ip2]) - d6*(by[i0]+by[ip1]);
            wp2 =  d3*(by[im1]+by[ip3]) + d5*(by[i0]+by[ip2]) - d4*by[ip1];
            wp3 =  -d2*(by[i0]+by[ip3]) + d1*(by[ip1]+by[ip2]);

            out[i0] = out[i0] + wm3*u[im3] + wm2*u[im2] + wm1*u[im1] + w0*u[i0] + wp1*u[ip1] + wp2*u[ip2] + wp3*u[ip3];
        }
    }
}


/* The gateway function */
/* u_x = D2x(u, bx, by, dx, dy, Nx, Ny, order, numThreads) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx, dy;          /* input scalar */
    double *u;          /* Nx1 input vector */
    double *bx, *by;          /* Nx1 variable coefficients */
    int Nx, Ny;         /* size of input */
    double *out;        /* output vector */
    int order;          /* Approximation order */

    /* create pointers to the input */
    u = mxGetPr(prhs[0]);
    bx = mxGetPr(prhs[1]);
    by = mxGetPr(prhs[2]);

    /* get scalars  */
    dx = mxGetScalar(prhs[3]);
    dy = mxGetScalar(prhs[4]);
    Nx = (int) mxGetScalar(prhs[5]);
    Ny = (int) mxGetScalar(prhs[6]);
    order = (int) mxGetScalar(prhs[7]);

    int numThreads = (int) mxGetScalar(prhs[8]);
    omp_set_num_threads(numThreads);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out = mxGetPr(plhs[0]);

    /* call the computational routine */
    if (order == 2){
        // D2_2nd(out, u, b, dx, Nx, Ny);
    }
    else if (order == 4){
        // D2_4th(out, u, b, dx, Nx, Ny);
    }
    else if (order == 6){
        D2_6th(out, u, bx, by, dx, dy, Nx, Ny);
    }
}