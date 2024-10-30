cpl		?= gfortran
f2py	?= f2py
equilibrium_solver ?= NLEQ1
PYTHON?=python3

SUPPRESS_WARNING_FOR_LEGACY_CODE = -w

# https://stackoverflow.com/questions/42419301/makefile-emit-warning-gcc-version-is-lower-than-4-8-0
GFORTRAN_LEGACYOPT = -std=legacy
GFORTRAN_LEGACY_VERSION = "8"
GFORTRAN_VERSION := "`gfortran -dumpversion`"
IS_GFORTRAN_ABOVE_LEGACY_VERSION := $(shell expr "$(GFORTRAN_VERSION)" ">=" "$(GFORTRAN_LEGACY_VERSION)")
ifeq ("$(IS_GFORTRAN_ABOVE_LEGACY_VERSION)", "0")
  GFORTRAN_LEGACYOPT =
endif

ifeq ($(cpl), ifort)
  lflags_debug = -debug -save-temps -fpic -heap-arrays -O0 -g -traceback -check all -fpe0 -fp-stack-check -fp-model precise -fimf-arch-consistency=true -fpp
  lflags_fast = -O2 -fpp
else
  lflags_debug = -ffpe-trap=invalid,zero,overflow -Og -fbacktrace -fcheck=all
  lflags_fast = -O2
endif

ifeq ($(equilibrium_solver), NLEQ2)
  eq_solver = nleq2
  eq_solver_switch = -DUSE_EQ_SOLVER_NLEQ2
else
  eq_solver = nleq1
  eq_solver_switch =
endif

LIB_DIR=.
MYRADEX_LIB=$(LIB_DIR)/libmy_radex.a

all: lflags = $(lflags_fast)
all: my_radex

wrapper: lflags = $(lflags_fast)
wrapper: wrapper_my_radex.so

debug: lflags = $(lflags_debug)
debug: my_radex

cflags = $(lflags) -c -cpp $(eq_solver_switch)

OBJSWRAPPER=linalg_nleq1.o my_radex.o nleq1.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o wnorm.o zibconst.o zibmon.o zibsec.o wrapper_for_python.o wrapper_for_cython.o pipe_fortran_python.o

cython_wrapper: myRadex.*.so

myRadex.*.so: setup.py wrapper_for_cython.pyx $(MYRADEX_LIB)
	$(PYTHON) setup.py build_ext --inplace

$(MYRADEX_LIB): $(OBJSWRAPPER)
	ar rcs $(MYRADEX_LIB) $(OBJSWRAPPER)

my_radex: configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o $(eq_solver).o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 
	$(cpl) $(lflags) -o my_radex configure.o main.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o $(eq_solver).o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 

wrapper_my_radex.so: wrapper_for_python.o my_radex.o opkda1.o opkda2.o opkdmain.o statistic_equilibrium.o sub_global_variables.o sub_trivials.o $(eq_solver).o linalg_nleq1.o wnorm.o zibconst.o zibmon.o zibsec.o 
	$(f2py) -c -m wrapper_my_radex wrapper_for_python.f90 --include-path "$(shell pwd)" --backend meson

configure.o: configure.f90 sub_trivials.o my_radex.o statistic_equilibrium.o
	$(cpl) $(cflags) configure.f90

main.o: main.f90 my_radex.o sub_trivials.o sub_global_variables.o statistic_equilibrium.o
	$(cpl) $(cflags) main.f90

my_radex.o: my_radex.f90 sub_trivials.o sub_global_variables.o statistic_equilibrium.o
	$(cpl) $(cflags) my_radex.f90

wrapper_for_cython.o: wrapper_for_cython.f90 wrapper_for_python.o sub_trivials.o sub_global_variables.o statistic_equilibrium.o
	$(cpl) $(cflags) wrapper_for_cython.f90

pipe_fortran_python.o: pipe_fortran_python.cpp pipe_fortran_python.hpp wrapper_for_cython.o
	g++ -std=c++11 -Wall -g -O3 -fPIC -c pipe_fortran_python.cpp

wrapper_for_python.o: wrapper_for_python.f90 my_radex.o
	$(cpl) $(cflags) wrapper_for_python.f90

opkda1.o:  opkda1.f
	$(cpl) $(cflags) $(GFORTRAN_LEGACYOPT) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) opkda1.f

opkda2.o:  opkda2.f
	$(cpl) $(cflags) $(GFORTRAN_LEGACYOPT) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) opkda2.f

opkdmain.o: opkdmain.f opkda1.o opkda2.o
	$(cpl) $(cflags) $(GFORTRAN_LEGACYOPT) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) opkdmain.f

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o $(eq_solver).o linalg_nleq1.o
	$(cpl) $(cflags) statistic_equilibrium.f90

$(eq_solver).o: $(eq_solver).f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) $(eq_solver).f

linalg_nleq1.o: linalg_nleq1.f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) linalg_nleq1.f

wnorm.o: wnorm.f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) wnorm.f

zibconst.o: zibconst.f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) zibconst.f

zibmon.o: zibmon.f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) zibmon.f

zibsec.o: zibsec.f
	$(cpl) $(cflags) $(SUPPRESS_WARNING_FOR_LEGACY_CODE) zibsec.f

sub_global_variables.o:  sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

sub_trivials.o: sub_trivials.f90 sub_global_variables.o
	$(cpl) $(cflags) sub_trivials.f90

clean:
	 rm ./*.mod ./*.o

CLEAN:
	 rm ./*.mod ./*.o ./*.so ./*.a
