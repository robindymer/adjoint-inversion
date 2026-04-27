
#include "mex.h"
#include <omp.h>

/* The computational routine */
// void D1_2nd(double *u_x, double *u, double dx, int Nx, int Ny)
// {
//     int i, j, bp;
//     double dxm = 1.0/(2*dx);
//     bp = 2;

//     for (j=bp; j<Nx-bp; j++){
//         for (i=0; i<Ny; i++) {
//             u_x[j*Ny + i] = dxm*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]);
//         }
//     }
// }

// void D1_4th(double *u_x, double *u, double dx, int Nx, int Ny)
// {
//     int i, j, bp;
//     double dxm = 1.0/dx;
//     double d[2] = {2.0/3.0, -1.0/12.0};
//     bp = 6;

//     for (j=bp; j<Nx-bp; j++){
//         for (i=0; i<Ny; i++) {
//             u_x[j*Ny + i] = dxm*(d[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + d[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]));
//         }
//     }
// }

void D1_6th(double *u_x, double*u_y, double *u, double dx, double dy, int Nx, int Ny)
{
    int i, j, bp;
    double dxm = 1.0/dx;
    double dym = 1.0/dy;
    double xd[3] = {dxm*3.0/4.0, -dxm*3.0/20.0, dxm*1.0/60.0};
    double yd[3] = {dym*3.0/4.0, -dym*3.0/20.0, dym*1.0/60.0};
    bp = 9;

    #pragma omp parallel for private(i,j)
    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            u_x[j*Ny + i] = xd[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + xd[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]) + xd[2]*(u[(j+3)*Ny + i] - u[(j-3)*Ny + i]);
            u_y[j*Ny + i] = yd[0]*(u[j*Ny + i+1] - u[j*Ny + i-1]) + yd[1]*(u[j*Ny + i+2] - u[j*Ny + i-2]) + yd[2]*(u[j*Ny + i+3] - u[j*Ny + i-3]);
        }
    }

    for (j=bp; j<Nx-bp; j++){
        for (i=0; i<bp; i++) {
            u_x[j*Ny + i] = xd[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + xd[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]) + xd[2]*(u[(j+3)*Ny + i] - u[(j-3)*Ny + i]);
        }
        for (i=Ny-bp; i<Ny; i++) {
            u_x[j*Ny + i] = xd[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + xd[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]) + xd[2]*(u[(j+3)*Ny + i] - u[(j-3)*Ny + i]);
        }
    }

    for (j=0; j<bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            u_y[j*Ny + i] = yd[0]*(u[j*Ny + i+1] - u[j*Ny + i-1]) + yd[1]*(u[j*Ny + i+2] - u[j*Ny + i-2]) + yd[2]*(u[j*Ny + i+3] - u[j*Ny + i-3]);
        }
    }

    for (j=Nx-bp; j<Nx; j++){
        for (i=bp; i<Ny-bp; i++) {
            u_y[j*Ny + i] = yd[0]*(u[j*Ny + i+1] - u[j*Ny + i-1]) + yd[1]*(u[j*Ny + i+2] - u[j*Ny + i-2]) + yd[2]*(u[j*Ny + i+3] - u[j*Ny + i-3]);
        }
    }

}

/* The gateway function */
/* [u_x, u_y] = D1x_and_D1y(u, dx, dy, Nx, Ny, order, numThreads) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx, dy;          /* input scalars */
    double *u;          /* Nx1 input matrix */
    int Nx, Ny;         /* size of matrix */
    double *out_x, *out_y;  /* output matrices */
    int order;          /* Approximation order */

    /* get grid spacing  */
    dx = mxGetScalar(prhs[1]);
    dy = mxGetScalar(prhs[2]);

    /* create a pointer to the input */
    u = mxGetPr(prhs[0]);

    /* get dimensions */
    Nx = (int) mxGetScalar(prhs[3]);
    Ny = (int) mxGetScalar(prhs[4]);
    order = (int) mxGetScalar(prhs[5]);

    int numThreads = (int) mxGetScalar(prhs[6]);
    omp_set_num_threads(numThreads);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out_x = mxGetPr(plhs[0]);
    out_y = mxGetPr(plhs[1]);

    /* call the computational routine */
    if (order == 2){
        // D1x_2nd(out, u, dx, Nx, Ny);
    }
    else if (order == 4){
        // D1x_4th(out, u, dx, Nx, Ny);
    }
    else if (order == 6){
        D1_6th(out_x, out_y, u, dx, dy, Nx, Ny);
    }
}