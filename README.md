**What's new**

- A function for calculating critical densities of all the transitions of a molecule is included.

---

This code solves essentially the same problem as
[```RADEX```](http://home.strw.leidenuniv.nl/~moldata/radex.html) written by
Van der Tak, except that we take a different approach to solve the statistical
equilibrium problem.  Given an initial distribution, what we do is to evolve
the system towards equilibrium using an ODE solver.

**Usage**

To use this code, you first need to compile it using the makefile (the
executable is named ```my_radex```) by running ```make``` in the command line,
then edit the configuration file (```configure.dat```) to meet your needs.
After this execute the command

    ./my_radex configure.dat

and you will get the results.

A python wrapper is also included.  To make the wrapper, run ```make
wrapper``` in the command line.  For its usage, see
[this Jupyter notebook](https://github.com/fjdu/myRadex/blob/master/example.ipynb).
The wrapper is preliminary; not all the functionalities in the Fortran source code are included in the wrapper
(though usually they are not needed).

We use the [LAMDA](http://home.strw.leidenuniv.nl/~moldata/molecules.html) format for the input energy levels and transition rates.

This code has not been thoroughly tested, and it is possible that some input files won't load properly.  Usually one can work around the problem by changing the input file to the format of a file that is known to work.

**Note**

- Works with python 3.7
