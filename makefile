#cpl	= mpif90
#cpl = ifort
#cflags = -debug -save-temps -fpic -Wl,-no_pie -heap-arrays -O2 -g -traceback -check all -fp-stack-check -c
#lflags = -debug -save-temps -fpic -heap-arrays -O2 -g -traceback -check all -fp-stack-check
#cflags = -debug -c
#lflags = -debug
#cflags = -fast -c
#lflags = -fast
cpl	= gfortran
lflags = -g -Wall -fcheck=all # -ffpe-trap=zero,invalid,overflow,underflow # -pedantic-errors
cflags = $(lflags) -c

OBJS = configure.o my_radex.o \
		main.o opkda1.o opkda2.o \
		opkdmain.o sub_global_variables.o sub_trivials.o \
		statistic_equilibrium.o

exe_name = a.out
exe_alt = a1.out

all: $(exe_name)

alt: $(exe_alt)

$(exe_name):  \
		$(OBJS)
	$(cpl) $(lflags) $(OBJS)

$(exe_alt):  \
		$(OBJS)
	$(cpl) $(lflags) -o $(exe_alt) $(OBJS)

main.o: main.f90 configure.o my_radex.o sub_trivials.o statistic_equilibrium.o
	$(cpl) $(cflags) main.f90

my_radex.o: my_radex.f90 sub_global_variables.o \
		statistic_equilibrium.o
	$(cpl) $(cflags) my_radex.f90

#data_struct.o: data_struct.f90
#	$(cpl) $(cflags) data_struct.f90

configure.o: configure.f90 my_radex.o sub_trivials.o
	$(cpl) $(cflags) configure.f90

sub_trivials.o: sub_trivials.f90
	$(cpl) $(cflags) sub_trivials.f90

sub_global_variables.o: sub_global_variables.f90
	$(cpl) $(cflags) sub_global_variables.f90

statistic_equilibrium.o: statistic_equilibrium.f90 sub_trivials.o sub_global_variables.o
	$(cpl) $(cflags) statistic_equilibrium.f90

opkdmain.o: opkdmain.f opkda1.o opkda2.o
	$(cpl) $(cflags) opkdmain.f

opkda1.o: opkda1.f
	$(cpl) $(cflags) opkda1.f

opkda2.o: opkda2.f
	$(cpl) $(cflags) opkda2.f

clean:
	 rm $(OBJS)
