To use this code, you first need to compile it using the makefile (the
executable is named a.out), then edit the configuration file (configure.dat) to
meet your needs.  After this execute the command
    ./a.out configure.dat.

A ```python``` wrapper is also included.  For its usage, see [this](http://nbviewer.ipython.org/urls/dl.dropbox.com/s/79pt2r35sqha51h/Myradex-Python-Wrapper-20140804.ipynb) ```ipython notebook```.  The wrapper is preliminary; not all the functionalities in the source code are included in the wrapper (though usually they are not needed).

We use the [LAMDA](http://home.strw.leidenuniv.nl/~moldata/molecules.html) format for the input energy levels and transition rates.

This code has not be thoroughly tested, and it is possible that some input file cannot be loaded properly.
