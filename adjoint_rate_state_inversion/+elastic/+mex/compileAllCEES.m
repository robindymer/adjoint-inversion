mex -outdir +elastic/+mex +elastic/+mex/D1x.c;
mex -outdir +elastic/+mex +elastic/+mex/D1y.c;
mex -outdir +elastic/+mex +elastic/+mex/D2xVariable.c;
mex -outdir +elastic/+mex +elastic/+mex/D2yVariable.c;
mex -outdir +elastic/+mex +elastic/+mex/D2xVariableDouble.c;
mex -outdir +elastic/+mex +elastic/+mex/D2yVariableDouble.c;

mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1xOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1yOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2xVariableOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2yVariableOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2xVariableDoubleOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2yVariableDoubleOMP.c;

mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1x_and_D1y_OMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2xVariablePlusD2yVariableOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1xAuPlusBvOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1yAuPlusBvOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D2VariableCombinedXPlusYOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1xAuPlusBvCombinedOMP.c;
mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1yAuPlusBvCombinedOMP.c;

mex CFLAGS="$CFLAGS -fopenmp -O3 -fwrapv" CLIBS="$CLIBS -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgomp" -outdir +elastic/+mex +elastic/+mex/D1LeftOMP.c;