To use this code, you first need to compile it using the makefile (the
executable is named ```my_radex```) by running ```make``` in the command line,
then edit the configuration file (```configure.dat```) to meet your needs.
After this execute the command
```./a.out configure.dat.```

A ```python``` wrapper is also included.  To make the wrapper, run ```make
wrapper``` in the command line.  For its usage, see
[this](http://nbviewer.ipython.org/urls/dl.dropbox.com/s/79pt2r35sqha51h/Myradex-Python-Wrapper-20140804.ipynb)
```ipython notebook```.  The wrapper is preliminary; not all the
functionalities in the ```Fortran``` source code are included in the wrapper
(though usually they are not needed).

We use the [LAMDA](http://home.strw.leidenuniv.nl/~moldata/molecules.html) format for the input energy levels and transition rates.

This code has not been thoroughly tested, and it is possible that some input files cannot be loaded properly.
