
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

    #pragma omp parallel for private(i,j)
    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            im3 = (j-3)*Ny + i;
            im2 = (j-2)*Ny + i;
            im1 = (j-1)*Ny + i;
            i0 = j*Ny + i;
            ip1 = (j+1)*Ny + i;
            ip2 = (j+2)*Ny + i;
            ip3 = (j+3)*Ny + i;

            out[i0] = d[0]*(bu[ip1]*u[ip1]-bu[im1]*u[im1]+bv[ip1]*v[ip1]-bv[im1]*v[im1]) + d[1]*(bu[ip2]*u[ip2]-bu[im2]*u[im2]+bv[ip2]*v[ip2]-bv[im2]*v[im2]) + d[2]*(bu[ip3]*u[ip3]-bu[im3]*u[im3]+bv[ip3]*v[ip3]-bv[im3]*v[im3]);
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