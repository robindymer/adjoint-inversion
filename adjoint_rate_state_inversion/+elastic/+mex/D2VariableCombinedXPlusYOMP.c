
#include "mex.h"
#include <omp.h>

/* The computational routines */
void D2x_2nd(double *u_x, double *v_x, double *u, double * v, double *b, double dx, int Nx, int Ny)
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
            v_x[idx] = dxm*(wm1*v[idx - str] + w0*v[idx] + wp1*v[idx + str]);
        }
    }
}

void D2x_4th(double *u_x, double *v_x, double *u, double *v, double *b, double dx, int Nx, int Ny)
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
            v_x[i0] = dxm*(wm2*v[im2] + wm1*v[im1] + w0*v[i0] + wp1*v[ip1] + wp2*v[ip2]);
        }
    }
}

void D2x_6th(double *u_x, double *v_x, double *u, double *v, double *bx, double *by, double dx, double dy, int Nx, int Ny)
{
    int i, j, bp, str, str2, str3;
    int im3, im2, im1, i0, ip1, ip2, ip3;
    double dxm = 1.0/(dx*dx);
    double dym = 1.0/(dy*dy);

    double wxm3, wxm2, wxm1, wx0, wxp1, wxp2, wxp3;
    double wym3, wym2, wym1, wy0, wyp1, wyp2, wyp3;
    double wyp3_m2, wyp3_m1, wyp2_m1;
    bp = 9;

    double cx1 = -dxm*1.0/40.0;
    double cx2 = -dxm*0.11e2/0.360e3;
    double cx3 = -dxm*1.0/0.20e2;
    double cx4 = -dxm*0.3e1/0.10e2;
    double cx5 = 0.7e1*cx1;
    double cx6 = 0.17e2*cx1;
    double cx7 = -dxm*1.0/0.180e3;
    double cx8 = -dxm*0.101e3/0.180e3;
    double cx9 = -dxm*1.0/0.8e1;
    double cx10 = 0.19e2*cx3;

    double cy1 = -dym*1.0/40.0;
    double cy2 = -dym*0.11e2/0.360e3;
    double cy3 = -dym*1.0/0.20e2;
    double cy4 = -dym*0.3e1/0.10e2;
    double cy5 = 0.7e1*cy1;
    double cy6 = 0.17e2*cy1;
    double cy7 = -dym*1.0/0.180e3;
    double cy8 = -dym*0.101e3/0.180e3;
    double cy9 = -dym*1.0/0.8e1;
    double cy10 = 0.19e2*cy3;

    str = Ny;
    str2 = 2*str;
    str3 = 3*str;

    #pragma omp parallel for private(i,j,im3,im2,im1,i0,ip1,ip2,ip3,wxm3,wxm2,wxm1,wx0,wxp1,wxp2,wxp3,wym3,wym2,wym1,wy0,wyp1,wyp2,wyp3,wyp3_m2,wyp3_m1,wyp2_m1)
    for (j=bp; j<Nx-bp; j++){

        i = bp;
        i0 = j*Ny + i;
        im3 = i0-3;
        im2 = i0-2;
        im1 = i0-1;
        ip1 = i0+1;
        ip2 = i0+2;
        ip3 = i0+3;

        wyp1 = -cy1*(by[im3]+by[ip2]) - cy4*(by[im2]+by[ip1]) - cy6*(by[im1]+by[i0]);

        wyp2_m1 = cy3*(by[im3]+by[ip1]) + cy5*(by[im2]+by[i0]) - cy4*by[im1];
        wyp2 = cy3*(by[im2]+by[ip2]) + cy5*(by[im1]+by[ip1]) - cy4*by[i0];

        wyp3_m2 = -cy2*(by[im3]+by[i0]) + cy1*(by[im2]+by[im1]);
        wyp3_m1 = -cy2*(by[im2]+by[ip1]) + cy1*(by[im1]+by[i0]);
        wyp3 = -cy2*(by[im1]+by[ip2]) + cy1*(by[i0]+by[ip1]);

        for (i=bp; i<Ny-bp; i++) {
            i0 = j*Ny + i;

            // Compute D2x
            im3 = i0-str3;
            im2 = i0-str2;
            im1 = i0-str;
            ip1 = i0+str;
            ip2 = i0+str2;
            ip3 = i0+str3;

            wxm3 = -cx2*(bx[im3]+bx[i0]) + cx1*(bx[im2]+bx[im1]);
            wxm2 = cx3*(bx[im3]+bx[ip1]) + cx5*(bx[im2]+bx[i0]) - cx4*bx[im1];
            wxm1 = -cx1*(bx[im3]+bx[ip2]) - cx4*(bx[im2]+bx[ip1]) - cx6*(bx[im1]+bx[i0]);
            wx0  =  cx7*(bx[im3]+bx[ip3]) + cx9*(bx[im2]+bx[ip2]) + cx10*(bx[im1]+bx[ip1]) + cx8*bx[i0];
            wxp1 = -cx1*(bx[im2]+bx[ip3]) - cx4*(bx[im1]+bx[ip2]) - cx6*(bx[i0]+bx[ip1]);
            wxp2 =  cx3*(bx[im1]+bx[ip3]) + cx5*(bx[i0]+bx[ip2]) - cx4*bx[ip1];
            wxp3 =  -cx2*(bx[i0]+bx[ip3]) + cx1*(bx[ip1]+bx[ip2]);

            u_x[i0] = wxm3*u[im3] + wxm2*u[im2] + wxm1*u[im1] + wx0*u[i0] + wxp1*u[ip1] + wxp2*u[ip2] + wxp3*u[ip3];
            v_x[i0] = wxm3*v[im3] + wxm2*v[im2] + wxm1*v[im1] + wx0*v[i0] + wxp1*v[ip1] + wxp2*v[ip2] + wxp3*v[ip3];

            // Compute D2y
            im3 = i0-3;
            im2 = i0-2;
            im1 = i0-1;
            ip1 = i0+1;
            ip2 = i0+2;
            ip3 = i0+3;

            wym3 = wyp3_m2;
            wyp3_m2 = wyp3_m1;
            wyp3_m1 = wyp3;

            wym2 = wyp2_m1;
            wyp2_m1 = wyp2;

            wym1 = wyp1;

            wy0 =   cy7*(by[im3]+by[ip3]) + cy9*(by[im2]+by[ip2]) + cy10*(by[im1]+by[ip1]) + cy8*by[i0];
            wyp1 = -cy1*(by[im2]+by[ip3]) - cy4*(by[im1]+by[ip2]) - cy6*(by[i0]+by[ip1]);
            wyp2 =  cy3*(by[im1]+by[ip3]) + cy5*(by[i0]+by[ip2]) - cy4*by[ip1];
            wyp3 =  -cy2*(by[i0]+by[ip3]) + cy1*(by[ip1]+by[ip2]);

            u_x[i0] += wym3*u[im3] + wym2*u[im2] + wym1*u[im1] + wy0*u[i0] + wyp1*u[ip1] + wyp2*u[ip2] + wyp3*u[ip3];
            v_x[i0] += wym3*v[im3] + wym2*v[im2] + wym1*v[im1] + wy0*v[i0] + wyp1*v[ip1] + wyp2*v[ip2] + wyp3*v[ip3];
        }
    }
}


/* The gateway function */
/* [u_x, v_x] = D2x(u, v, b, dx, Nx, Ny, order) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double dx, dy;          /* input scalar */
    double *u, *v;      /* Nx1 input vectors */
    double *bx, *by;          /* Nx1 variable coefficient */
    int Nx, Ny;         /* size of input */
    double *out1, *out2;        /* output vector */
    int order;          /* Approximation order */

    /* create pointers to the input */
    u = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    bx = mxGetPr(prhs[2]);
    by = mxGetPr(prhs[3]);

    /* get scalars  */
    dx = mxGetScalar(prhs[4]);
    dy = mxGetScalar(prhs[5]);
    Nx = (int) mxGetScalar(prhs[6]);
    Ny = (int) mxGetScalar(prhs[7]);
    order = (int) mxGetScalar(prhs[8]);

    int numThreads = (int) mxGetScalar(prhs[9]);
    omp_set_num_threads(numThreads);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Nx*Ny, 1, mxREAL);

    /* get a pointer to the real data in the output matrix */
    out1 = mxGetPr(plhs[0]);
    out2 = mxGetPr(plhs[1]);

    /* call the computational routine */
    if (order == 2){
        // D2x_2nd(out1, out2, u, v, b, dx, Nx, Ny);
    }
    else if (order == 4){
        // D2x_4th(out1, out2, u, v, b, dx, Nx, Ny);
    }
    else if (order == 6){
        D2x_6th(out1, out2, u, v, bx, by, dx, dy, Nx, Ny);
    }
}