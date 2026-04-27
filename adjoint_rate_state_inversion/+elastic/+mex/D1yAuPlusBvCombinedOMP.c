
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

void D1x_6th(double *out1, double *out2, double *u, double *v, double *bu1, double *bv1, double *bu2, double *bv2, double dx, int Nx, int Ny)
{
    int i, j, bp, im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/dx;
    double d[3] = {dxm*3.0/4.0, -dxm*3.0/20.0, dxm*1.0/60.0};
    bp = 9;

    double x1_im3, x1_im2, x1_im1, x1_i0, x1_ip1, x1_ip2, x1_ip3;
    double x2_im3, x2_im2, x2_im1, x2_i0, x2_ip1, x2_ip2, x2_ip3;

    #pragma omp parallel for private(i,j, im3, im2, im1, i0, ip1, ip2, ip3, x1_im3, x1_im2, x1_im1, x1_i0, x1_ip1, x1_ip2, x1_ip3, x2_im3, x2_im2, x2_im1, x2_i0, x2_ip1, x2_ip2, x2_ip3)
    for (j=bp; j<Nx-bp; j++){
        i = bp-1;

        im2 = j*Ny + i-2;
        im1 = j*Ny + i-1;
        i0  = j*Ny + i;
        ip1 = j*Ny + i+1;
        ip2 = j*Ny + i+2;
        ip3 = j*Ny + i+3;

        x1_im2 = bu1[im2]*u[im2] + bv1[im2]*v[im2];
        x1_im1 = bu1[im1]*u[im1] + bv1[im1]*v[im1];
        x1_i0 = bu1[i0]*u[i0] + bv1[i0]*v[i0];
        x1_ip1 = bu1[ip1]*u[ip1] + bv1[ip1]*v[ip1];
        x1_ip2 = bu1[ip2]*u[ip2] + bv1[ip2]*v[ip2];
        x1_ip3 = bu1[ip3]*u[ip3] + bv1[ip3]*v[ip3];

        x2_im2 = bu2[im2]*u[im2] + bv2[im2]*v[im2];
        x2_im1 = bu2[im1]*u[im1] + bv2[im1]*v[im1];
        x2_i0 = bu2[i0]*u[i0] + bv2[i0]*v[i0];
        x2_ip1 = bu2[ip1]*u[ip1] + bv2[ip1]*v[ip1];
        x2_ip2 = bu2[ip2]*u[ip2] + bv2[ip2]*v[ip2];
        x2_ip3 = bu2[ip3]*u[ip3] + bv2[ip3]*v[ip3];

        for (i=bp; i<Ny-bp; i++) {

            // First output
            x1_im3 = x1_im2;
            x1_im2 = x1_im1;
            x1_im1 = x1_i0;
            x1_i0 = x1_ip1;
            x1_ip1 = x1_ip2;
            x1_ip2 = x1_ip3;

            ip3 = j*Ny + i+3;
            x1_ip3 = bu1[ip3]*u[ip3] + bv1[ip3]*v[ip3];

            i0  = j*Ny + i;
            out1[i0] = d[0]*(x1_ip1 - x1_im1) + d[1]*(x1_ip2 - x1_im2) + d[2]*(x1_ip3 - x1_im3);

            // Second output
            x2_im3 = x2_im2;
            x2_im2 = x2_im1;
            x2_im1 = x2_i0;
            x2_i0 = x2_ip1;
            x2_ip1 = x2_ip2;
            x2_ip2 = x2_ip3;

            x2_ip3 = bu2[ip3]*u[ip3] + bv2[ip3]*v[ip3];

            out2[i0] = d[0]*(x2_ip1 - x2_im1) + d[1]*(x2_ip2 - x2_im2) + d[2]*(x2_ip3 - x2_im3);
        }
    }
}

/* The gateway function */
/* out = D1yAuPlusBvCombinedOMP(u, v, bu1, bv1, bu2, bv2, dx, Nx, Ny, order, numThreads) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx;          /* input scalar */
    double *u, *v, *bu1, *bv1, *bu2, *bv2;  /* Nx1 input vectors */
    int Nx, Ny;         /* size of matrix */
    double *out1, *out2;        /* output vectors */
    int order;          /* Approximation order */



    /* create pointers to the input */
    u = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    bu1 = mxGetPr(prhs[2]);
    bv1 = mxGetPr(prhs[3]);
    bu2 = mxGetPr(prhs[4]);
    bv2 = mxGetPr(prhs[5]);

    /* get grid spacing  */
    dx = mxGetScalar(prhs[6]);

    /* get dimensions */
    Nx = (int) mxGetScalar(prhs[7]);
    Ny = (int) mxGetScalar(prhs[8]);
    order = (int) mxGetScalar(prhs[9]);

    int numThreads = (int) mxGetScalar(prhs[10]);
    omp_set_num_threads(numThreads);

    /* create the output vectors */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out1 = mxGetPr(plhs[0]);
    out2 = mxGetPr(plhs[1]);

    /* call the computational routine */
    if (order == 2){
        // D1x_2nd(out, u, dx, Nx, Ny);
    }
    else if (order == 4){
        // D1x_4th(out, u, dx, Nx, Ny);
    }
    else if (order == 6){
        D1x_6th(out1, out2, u, v, bu1, bv1, bu2, bv2, dx, Nx, Ny);
    }
}