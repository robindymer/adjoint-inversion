
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

void D1x_6th(double *out1, double *out2, double *ux, double *vx, double *uy, double *vy, double *bu1, double *bv1, double *bu2, double *bv2, double dx, double dy, int Nx, int Ny)
{
    int i, j, bp, im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/dx;
    double dym = 1.0/dy;
    double xd[3] = {dxm*3.0/4.0, -dxm*3.0/20.0, dxm*1.0/60.0};
    double yd[3] = {dym*3.0/4.0, -dym*3.0/20.0, dym*1.0/60.0};
    bp = 9;

    double x1_im3, x1_im2, x1_im1, x1_i0, x1_ip1, x1_ip2, x1_ip3;
    double x2_im3, x2_im2, x2_im1, x2_i0, x2_ip1, x2_ip2, x2_ip3;

    #pragma omp parallel for private(i,j, im3, im2, im1, i0, ip1, ip2, ip3, x1_im3, x1_im2, x1_im1, x1_i0, x1_ip1, x1_ip2, x1_ip3, x2_im3, x2_im2, x2_im1, x2_i0, x2_ip1, x2_ip2, x2_ip3)
    for (j=bp; j<Nx-bp; j++){

        // Setup Dy
        i = bp-1;

        im2 = j*Ny + i-2;
        im1 = j*Ny + i-1;
        i0  = j*Ny + i;
        ip1 = j*Ny + i+1;
        ip2 = j*Ny + i+2;
        ip3 = j*Ny + i+3;

        x1_im2 = bu1[im2]*uy[im2] + bu2[im2]*vy[im2];
        x1_im1 = bu1[im1]*uy[im1] + bu2[im1]*vy[im1];
        x1_i0 = bu1[i0]*uy[i0] + bu2[i0]*vy[i0];
        x1_ip1 = bu1[ip1]*uy[ip1] + bu2[ip1]*vy[ip1];
        x1_ip2 = bu1[ip2]*uy[ip2] + bu2[ip2]*vy[ip2];
        x1_ip3 = bu1[ip3]*uy[ip3] + bu2[ip3]*vy[ip3];

        x2_im2 = bv1[im2]*uy[im2] + bv2[im2]*vy[im2];
        x2_im1 = bv1[im1]*uy[im1] + bv2[im1]*vy[im1];
        x2_i0 = bv1[i0]*uy[i0] + bv2[i0]*vy[i0];
        x2_ip1 = bv1[ip1]*uy[ip1] + bv2[ip1]*vy[ip1];
        x2_ip2 = bv1[ip2]*uy[ip2] + bv2[ip2]*vy[ip2];
        x2_ip3 = bv1[ip3]*uy[ip3] + bv2[ip3]*vy[ip3];

        for (i=bp; i<Ny-bp; i++) {

            // Dx
            im3 = (j-3)*Ny + i;
            im2 = (j-2)*Ny + i;
            im1 = (j-1)*Ny + i;
            i0 = j*Ny + i;
            ip1 = (j+1)*Ny + i;
            ip2 = (j+2)*Ny + i;
            ip3 = (j+3)*Ny + i;

            out1[i0] = xd[0]*(bu1[ip1]*ux[ip1]-bu1[im1]*ux[im1]+bv1[ip1]*vx[ip1]-bv1[im1]*vx[im1]) + xd[1]*(bu1[ip2]*ux[ip2]-bu1[im2]*ux[im2]+bv1[ip2]*vx[ip2]-bv1[im2]*vx[im2]) + xd[2]*(bu1[ip3]*ux[ip3]-bu1[im3]*ux[im3]+bv1[ip3]*vx[ip3]-bv1[im3]*vx[im3]);
            out2[i0] = xd[0]*(bu2[ip1]*ux[ip1]-bu2[im1]*ux[im1]+bv2[ip1]*vx[ip1]-bv2[im1]*vx[im1]) + xd[1]*(bu2[ip2]*ux[ip2]-bu2[im2]*ux[im2]+bv2[ip2]*vx[ip2]-bv2[im2]*vx[im2]) + xd[2]*(bu2[ip3]*ux[ip3]-bu2[im3]*ux[im3]+bv2[ip3]*vx[ip3]-bv2[im3]*vx[im3]);

            // Dy
            // First output
            x1_im3 = x1_im2;
            x1_im2 = x1_im1;
            x1_im1 = x1_i0;
            x1_i0 = x1_ip1;
            x1_ip1 = x1_ip2;
            x1_ip2 = x1_ip3;

            ip3 = j*Ny + i+3;
            x1_ip3 = bu1[ip3]*uy[ip3] + bu2[ip3]*vy[ip3];

            out1[i0] += yd[0]*(x1_ip1 - x1_im1) + yd[1]*(x1_ip2 - x1_im2) + yd[2]*(x1_ip3 - x1_im3);

            // Second output
            x2_im3 = x2_im2;
            x2_im2 = x2_im1;
            x2_im1 = x2_i0;
            x2_i0 = x2_ip1;
            x2_ip1 = x2_ip2;
            x2_ip2 = x2_ip3;

            x2_ip3 = bv1[ip3]*uy[ip3] + bv2[ip3]*vy[ip3];

            out2[i0] += yd[0]*(x2_ip1 - x2_im1) + yd[1]*(x2_ip2 - x2_im2) + yd[2]*(x2_ip3 - x2_im3);
        }
    }
}

/* The gateway function */
/* out = D1Left(ux, vx, uy, vy, bu1, bv1, bu2, bv2, dx, Nx, Ny, order, numThreads) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx, dy;          /* input scalar */
    double *ux, *vx, *uy, *vy, *bu1, *bv1, *bu2, *bv2;  /* Nx1 input vectors */
    int Nx, Ny;         /* size of matrix */
    double *out1, *out2;        /* output vectors */
    int order;          /* Approximation order */

    /* create pointers to the input */
    ux = mxGetPr(prhs[0]);
    vx = mxGetPr(prhs[1]);
    uy = mxGetPr(prhs[2]);
    vy = mxGetPr(prhs[3]);
    bu1 = mxGetPr(prhs[4]);
    bv1 = mxGetPr(prhs[5]);
    bu2 = mxGetPr(prhs[6]);
    bv2 = mxGetPr(prhs[7]);

    /* get grid spacing  */
    dx = mxGetScalar(prhs[8]);
    dy = mxGetScalar(prhs[9]);

    /* get dimensions */
    Nx = (int) mxGetScalar(prhs[10]);
    Ny = (int) mxGetScalar(prhs[11]);
    order = (int) mxGetScalar(prhs[12]);

    int numThreads = (int) mxGetScalar(prhs[13]);
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
        D1x_6th(out1, out2, ux, vx, uy, vy, bu1, bv1, bu2, bv2, dx, dy, Nx, Ny);
    }
}