cpl	?= gfortran

all: my_radex

wrapper: wrapper_my_radex.so


lflags = -O3 -fPIC
cflags = $(lflags) -c

my_radex: configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o
	$(cpl) $(lflags) -o my_radex configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o

wrapper_my_radex.so: wrapper_for_python.f90 my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o
	f2py -c -m wrapper_my_radex wrapper_for_python.f90 -L my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o

configure.o: configure.f90 sub_trivials.o my_radex.o statistic_equilibrium.o
	$(cpl) $(cflags) configure.f90

main.o: main.f90 my_radex.o sub_trivials.o sub_global_variables.o statistic_equilibrium.o
	$(cpl) $(cflags) main.f90

my_radex.o: my_radex.f90 sub_trivials.o sub_global_variables.o statistic_equilibrium.o
	$(cpl) $(cflags) my_radex.f90

opkda1.o:  opkda1.f
	$(cpl) $(cflags) opkda1.f

opkda2.o:  opkda2.f
	$(cpl) $(cflags) opkda2.f

opkdmain.o: opkdmain.f opkda1.o opkda2.o
	$(cpl) $(cflags) opkdmain.f

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o
	$(cpl) $(cflags) statistic_equilibrium.f90

sub_global_variables.o:  sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

sub_trivials.o:          sub_trivials.f90 sub_global_variables.o
	$(cpl) $(cflags) sub_trivials.f90

clean:
	 rm ./*.mod ./*.o
