**What's new**

- 2023-03-29: The critical densities (calculated in two manners) are included in the output file.  Note that different collisional partners have different sets of critical densities.  At present only the values for the first collisional partner are included in the output file.

- 2022-06-06: `myRadex` is now included in the Astrophysics Source Code Library as ["myRadex: Radex with a twist"](https://ascl.net/2205.011).

- 2021-11-25: Instead of being calculated from the energy levels, now the frequencies in the input file will be used by default.  Also the energy level numbers will not be subtracted by 1 by default.  For backward compatibility, two Boolean options are added: `recalculateFreqWithEupElow` and `iLevel_subtract_one` (both are `False` by default).
- 2020-04-05: A function for calculating critical densities of all the transitions of a molecule is included.

---

## Calculation of critical density

By default equation (9) of [Shirley2015](https://iopscience.iop.org/article/10.1086/680342) is used.

A simpler definition that is only valid for two-level system is also used (labeled as "n_crit_old" in the output).  This is described in the text below equation (5) of [Shirley2015](https://iopscience.iop.org/article/10.1086/680342).

---

This code solves essentially the same problem as
[`RADEX`](http://home.strw.leidenuniv.nl/~moldata/radex.html) written by Van der Tak, except that we take a different approach to solve the statistical equilibrium problem.
Given an initial distribution, what we do is to evolve the system towards equilibrium using an ODE solver.

## Usage

To use this code, you first need to compile it using the `makefile` (the executable is named `my_radex`) by running `make` in the command line, then edit the configuration file (`configure.dat`) to meet your needs.  After this execute the command
```
./my_radex configure.dat
```
and you will get the results.

A python wrapper is also included.  To make the wrapper, run `make wrapper` in the command line.
For its usage, see [this Jupyter notebook](https://github.com/fjdu/myRadex/blob/master/example.ipynb).
The wrapper is preliminary; not all the functionalities in the Fortran source code are included in the wrapper (though usually they are not needed).

We use the [LAMDA](http://home.strw.leidenuniv.nl/~moldata/molecules.html) format for the input energy levels and transition rates.

This code has not been thoroughly tested, and it is possible that some input files won't load properly.
Usually one can work around the problem by changing the input file to the format of a file that is known to work.
