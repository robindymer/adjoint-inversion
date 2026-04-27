
#include "mex.h"
#include <omp.h>

/* The computational routine */
void D1x_2nd(double *u_x, double *u, double dx, int Nx, int Ny)
{
    int i, j, bp;
    double dxm = 1.0/(2*dx);
    bp = 2;

    for (j=bp; j<Nx-bp; j++){
        for (i=0; i<Ny; i++) {
            u_x[j*Ny + i] = dxm*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]);
        }
    }
}

void D1x_4th(double *u_x, double *u, double dx, int Nx, int Ny)
{
    int i, j, bp;
    double dxm = 1.0/dx;
    double d[2] = {2.0/3.0, -1.0/12.0};
    bp = 6;

    for (j=bp; j<Nx-bp; j++){
        for (i=0; i<Ny; i++) {
            u_x[j*Ny + i] = dxm*(d[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + d[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]));
        }
    }
}

void D1x_6th(double *out, double *u, double *v, double *bu, double *bv, double dx, int Nx, int Ny)
{
    int i, j, bp, im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/dx;
    double d[3] = {dxm*3.0/4.0, -dxm*3.0/20.0, dxm*1.0/60.0};
    bp = 9;

    double x_im3, x_im2, x_im1, x_i0, x_ip1, x_ip2, x_ip3;

    #pragma omp parallel for private(i,j, im3, im2, im1, i0, ip1, ip2, ip3, x_im3, x_im2, x_im1, x_i0, x_ip1, x_ip2, x_ip3)
    for (j=bp; j<Nx-bp; j++){
        i = bp-1;

        im2 = j*Ny + i-2;
        im1 = j*Ny + i-1;
        i0  = j*Ny + i;
        ip1 = j*Ny + i+1;
        ip2 = j*Ny + i+2;
        ip3 = j*Ny + i+3;

        x_im2 = bu[im2]*u[im2] + bv[im2]*v[im2];
        x_im1 = bu[im1]*u[im1] + bv[im1]*v[im1];
        x_i0 = bu[i0]*u[i0] + bv[i0]*v[i0];
        x_ip1 = bu[ip1]*u[ip1] + bv[ip1]*v[ip1];
        x_ip2 = bu[ip2]*u[ip2] + bv[ip2]*v[ip2];
        x_ip3 = bu[ip3]*u[ip3] + bv[ip3]*v[ip3];

        for (i=bp; i<Ny-bp; i++) {
            x_im3 = x_im2;
            x_im2 = x_im1;
            x_im1 = x_i0;
            x_i0 = x_ip1;
            x_ip1 = x_ip2;
            x_ip2 = x_ip3;

            ip3 = j*Ny + i+3;
            x_ip3 = bu[ip3]*u[ip3] + bv[ip3]*v[ip3];

            i0  = j*Ny + i;
            out[i0] = d[0]*(x_ip1 - x_im1) + d[1]*(x_ip2 - x_im2) + d[2]*(x_ip3 - x_im3);
        }
    }
}

/* The gateway function */
/* out = D1xAuPlusBvOMP(u, v, bu, bv, dx, Nx, Ny, order, numThreads) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx;          /* input scalar */
    double *u, *v, *bu, *bv;  /* Nx1 input vectors */
    int Nx, Ny;         /* size of matrix */
    double *out;        /* output matrix */
    int order;          /* Approximation order */



    /* create pointers to the input */
    u = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    bu = mxGetPr(prhs[2]);
    bv = mxGetPr(prhs[3]);

    /* get grid spacing  */
    dx = mxGetScalar(prhs[4]);

    /* get dimensions */
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
        // D1x_2nd(out, u, dx, Nx, Ny);
    }
    else if (order == 4){
        // D1x_4th(out, u, dx, Nx, Ny);
    }
    else if (order == 6){
        D1x_6th(out, u, v, bu, bv, dx, Nx, Ny);
    }
}