all:
	/usr/local/bin/f2py -c -m gadgetlibompv2 gadgetlibv2.f90 --f90flags="-fopenmp"  -lgomp
	#

clean:
	rm gadgetlibompv2.o gadgetlibompv2.so
