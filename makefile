cpl		?= gfortran
f2py	?= f2py

ifeq ($(cpl), ifort)
  lflags_debug = -debug -save-temps -fpic -heap-arrays -O0 -g -traceback -check all -fpe0 -fp-stack-check -fp-model precise -fimf-arch-consistency=true -fpp
  lflags_fast = -O2 -fpp
else
  lflags_debug = -ffpe-trap=invalid,zero,overflow -Og -fbacktrace -fcheck=all
  lflags_fast = -O2
endif

all: lflags = $(lflags_fast)
all: my_radex

wrapper: lflags = $(lflags_fast)
wrapper: wrapper_my_radex.so

debug: lflags = $(lflags_debug)
debug: my_radex

cflags = $(lflags) -c

my_radex: configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o nleq1.o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 
	$(cpl) $(lflags) -o my_radex configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o nleq1.o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 

wrapper_my_radex.so: wrapper_for_python.f90 my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o nleq1.o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 
	$(f2py) -c -m wrapper_my_radex wrapper_for_python.f90 my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o nleq1.o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 

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

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o nleq1.o linalg_nleq1.o
	$(cpl) $(cflags) statistic_equilibrium.f90

nleq1.o: nleq1.f
	$(cpl) $(cflags) nleq1.f

linalg_nleq1.o: linalg_nleq1.f
	$(cpl) $(cflags) linalg_nleq1.f

wnorm.o: wnorm.f
	$(cpl) $(cflags) wnorm.f

zibconst.o: zibconst.f
	$(cpl) $(cflags) zibconst.f

zibmon.o: zibmon.f
	$(cpl) $(cflags) zibmon.f

zibsec.o: zibsec.f
	$(cpl) $(cflags) zibsec.f

sub_global_variables.o:  sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

sub_trivials.o: sub_trivials.f90 sub_global_variables.o
	$(cpl) $(cflags) sub_trivials.f90

clean:
	 rm ./*.mod ./*.o

CLEAN:
	 rm ./*.mod ./*.o ./*.so ./*.a
