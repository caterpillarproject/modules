rm gadgetlibomp.o gadgetlibomp.so
f2py -c -m gadgetlibomp gadgetlib.f90 --fcompiler=gfortran --f90flags="-fopenmp" -lgomp
ipython