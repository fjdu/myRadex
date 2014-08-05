all: wrapper_my_radex.so

wrapper_my_radex.so: wrapper_for_python.f90 statistic_equilibrium.f90 my_radex.f90 opkda1.f opkda2.f opkdmain.f
	f2py-2.7 -c -m wrapper_my_radex wrapper_for_python.f90 -I"./" -L *.o # statistic_equilibrium.o my_radex.o opkda1.o opkda2.o opkdmain.o
