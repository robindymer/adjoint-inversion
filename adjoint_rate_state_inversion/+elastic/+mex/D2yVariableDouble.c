
#include "mex.h"

/* The computational routines */
void D2y_2nd(double *u_x, double *v_x, double *u, double *v, double *b, double dx, int Nx, int Ny)
{
    int i, j, bp, idx;
    double dxm = 1.0/(dx*dx);

    double wm1, w0, wp1;
    bp = 2;

    for (j=bp; j<Nx-bp; j++){
        for (i=bp; i<Ny-bp; i++) {
            idx = j*Ny + i;

            wm1 = 0.5*(b[idx - 1] + b[idx]);
            w0 = -0.5*(b[idx - 1] + b[idx + 1]) - b[idx];
            wp1 = 0.5*(b[idx + 1] + b[idx]);

            u_x[idx] = dxm*(wm1*u[idx - 1] + w0*u[idx] + wp1*u[idx + 1]);
            v_x[idx] = dxm*(wm1*v[idx - 1] + w0*v[idx] + wp1*v[idx + 1]);
        }
    }
}

void D2y_4th(double *u_x, double *v_x, double *u, double *v, double *b, double dx, int Nx, int Ny)
{
    int i, j, bp;
    int im2, im1, i0, ip1, ip2;
    double dxm = 1.0/(dx*dx);

    double wm2, wm1, w0, wp1, wp2, wp2_old;
    bp = 6;

    double c1 = 1.0/8.0;
    double c2 = 1.0/6.0;
    double c3 = 1.0/24.0;
    double c4 = 1.0/1.2;

    for (j=bp; j<Nx-bp; j++){
        i = bp;
        i0 = j*Ny + i;
        im2 = i0-2;
        im1 = i0-1;
        ip1 = i0+1;
        ip2 = i0+2;

        wp1 = (b[im2] + b[ip1])*c2 + (b[im1] + b[i0])*0.5;

        wp2_old = -(b[im2] + b[i0])*c1 + b[im1]*c2;
        wp2 = -(b[im1] + b[ip1])*c1 + b[i0]*c2;

        for (i=bp; i<Ny-bp; i++) {
            i0 = j*Ny + i;
            im2 = i0-2;
            im1 = i0-1;
            ip1 = i0+1;
            ip2 = i0+2;

            // wm2 = -(b[im2] + b[i0])*c1 + b[im1]*c2;
            wm2 = wp2_old;
            wp2_old = wp2;
            // wm1 = (b[im2] + b[ip1])*c2 + (b[im1] + b[i0])*0.5;
            wm1 = wp1;
            w0 = -(b[im2] + b[ip2])*c3 -(b[im1] + b[ip1])*c4 - 0.75*b[i0];
            wp1 = (b[im1] + b[ip2])*c2 + (b[i0] + b[ip1])*0.5;
            wp2 = -(b[i0] + b[ip2])*c1 + b[ip1]*c2;

            u_x[i0] = dxm*(wm2*u[im2] + wm1*u[im1] + w0*u[i0] + wp1*u[ip1] + wp2*u[ip2]);
            v_x[i0] = dxm*(wm2*v[im2] + wm1*v[im1] + w0*v[i0] + wp1*v[ip1] + wp2*v[ip2]);
        }
    }
}

void D2y_6th(double *u_x, double *v_x, double *u, double *v, double *b, double dx, int Nx, int Ny)
{
    int i, j, bp;
    int im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/(dx*dx);

    double wm3, wm2, wm1, w0, wp1, wp2, wp3;
    double wp2_m1, wp3_m1, wp3_m2;
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

    for (j=bp; j<Nx-bp; j++){
        i = bp;
        i0 = j*Ny + i;
        im3 = i0-3;
        im2 = i0-2;
        im1 = i0-1;
        ip1 = i0+1;
        ip2 = i0+2;
        ip3 = i0+3;

        wp1 = -c1*(b[im3]+b[ip2]) - c4*(b[im2]+b[ip1]) - c6*(b[im1]+b[i0]);

        wp2_m1 = c3*(b[im3]+b[ip1]) + c5*(b[im2]+b[i0]) - c4*b[im1];
        wp2 = c3*(b[im2]+b[ip2]) + c5*(b[im1]+b[ip1]) - c4*b[i0];

        wp3_m2 = -c2*(b[im3]+b[i0]) + c1*(b[im2]+b[im1]);
        wp3_m1 = -c2*(b[im2]+b[ip1]) + c1*(b[im1]+b[i0]);
        wp3 = -c2*(b[im1]+b[ip2]) + c1*(b[i0]+b[ip1]);

        for (i=bp; i<Ny-bp; i++) {
            i0 = j*Ny + i;
            im3 = i0-3;
            im2 = i0-2;
            im1 = i0-1;
            ip1 = i0+1;
            ip2 = i0+2;
            ip3 = i0+3;

            // wm3 = -c2*(b[im3]+b[i0]) + c1*(b[im2]+b[im1]);
            wm3 = wp3_m2;
            wp3_m2 = wp3_m1;
            wp3_m1 = wp3;

            // wm2 = c3*(b[im3]+b[ip1]) + c5*(b[im2]+b[i0]) - c4*b[im1];
            wm2 = wp2_m1;
            wp2_m1 = wp2;

            // wm1 = -c1*(b[im3]+b[ip2]) - c4*(b[im2]+b[ip1]) - c6*(b[im1]+b[i0]);
            wm1 = wp1;

            w0 =   c7*(b[im3]+b[ip3]) + c9*(b[im2]+b[ip2]) + c10*(b[im1]+b[ip1]) + c8*b[i0];
            wp1 = -c1*(b[im2]+b[ip3]) - c4*(b[im1]+b[ip2]) - c6*(b[i0]+b[ip1]);
            wp2 =  c3*(b[im1]+b[ip3]) + c5*(b[i0]+b[ip2]) - c4*b[ip1];
            wp3 =  -c2*(b[i0]+b[ip3]) + c1*(b[ip1]+b[ip2]);

            u_x[i0] = wm3*u[im3] + wm2*u[im2] + wm1*u[im1] + w0*u[i0] + wp1*u[ip1] + wp2*u[ip2] + wp3*u[ip3];
            v_x[i0] = wm3*v[im3] + wm2*v[im2] + wm1*v[im1] + w0*v[i0] + wp1*v[ip1] + wp2*v[ip2] + wp3*v[ip3];
        }
    }
}


/* The gateway function */
/* u_x = D2y(u, v, b, dx, Nx, Ny, order) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx;          /* input scalar */
    double *u, *v;      /* Nx1 input vectors */
    double *b;          /* Nx1 variable coefficient */
    int Nx, Ny;         /* size of input */
    double *out1, *out2; /* output vectors */
    int order;          /* Approximation order */

    /* create pointers to the input */
    u = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    b = mxGetPr(prhs[2]);

    /* get scalars  */
    dx = mxGetScalar(prhs[3]);
    Nx = (int) mxGetScalar(prhs[4]);
    Ny = (int) mxGetScalar(prhs[5]);
    order = (int) mxGetScalar(prhs[6]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out1 = mxGetPr(plhs[0]);
    out2 = mxGetPr(plhs[1]);

    /* call the computational routine */
    if (order == 2){
        D2y_2nd(out1, out2, u, v, b, dx, Nx, Ny);
    }
    else if (order == 4){
        D2y_4th(out1, out2, u, v, b, dx, Nx, Ny);
    }
    else if (order == 6){
        D2y_6th(out1, out2, u, v, b, dx, Nx, Ny);
    }
}