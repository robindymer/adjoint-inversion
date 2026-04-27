
#include "mex.h"

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

void D1x_6th(double *u_x, double *u, double dx, int Nx, int Ny)
{
    int i, j, bp;
    double dxm = 1.0/dx;
    double d[3] = {dxm*3.0/4.0, -dxm*3.0/20.0, dxm*1.0/60.0};
    bp = 9;

    for (j=bp; j<Nx-bp; j++){
        for (i=0; i<Ny; i++) {
            u_x[j*Ny + i] = d[0]*(u[(j+1)*Ny + i] - u[(j-1)*Ny + i]) + d[1]*(u[(j+2)*Ny + i] - u[(j-2)*Ny + i]) + d[2]*(u[(j+3)*Ny + i] - u[(j-3)*Ny + i]);
        }
    }
}

/* The gateway function */
/* u_x = D1x(u, dx, Nx, Ny, order) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx;          /* input scalar */
    double *u;          /* Nx1 input matrix */
    int Nx, Ny;         /* size of matrix */
    double *out;        /* output matrix */
    int order;          /* Approximation order */

    /* get grid spacing  */
    dx = mxGetScalar(prhs[1]);

    /* create a pointer to the input */
    u = mxGetPr(prhs[0]);

    /* get dimensions */
    Nx = (int) mxGetScalar(prhs[2]);
    Ny = (int) mxGetScalar(prhs[3]);
    order = (int) mxGetScalar(prhs[4]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out = mxGetPr(plhs[0]);

    /* call the computational routine */
    if (order == 2){
        D1x_2nd(out, u, dx, Nx, Ny);
    }
    else if (order == 4){
        D1x_4th(out, u, dx, Nx, Ny);
    }
    else if (order == 6){
        D1x_6th(out, u, dx, Nx, Ny);
    }
}