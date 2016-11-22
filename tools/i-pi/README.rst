i-PI V1.0 -- LAMMPS 
-------------------

A Python interface for ab initio path integral molecular dynamics simulations. 
i-PI is composed of a Python server (i-pi itself, that does not need to be 
compiled but only requires a relatively recent version of Python and Numpy)
that propagates the (path integral) dynamics of the nuclei, and of an external
code that acts as client and computes the electronic energy and forces.

This is typically a patched version of an electronic structure code, but a 
simple self-contained Fortran driver that implements Lennard-Jones and 
Silveira-Goldman potentials is included for test purposes.

This folder contains a stripped-down version to be used with LAMMPS, and might
not contain all the features of the latest version. Please see 
[http://epfl-cosmo.github.io/gle4md/index.html?page=ipi] or
[http://github.com/i-pi/i-pi] to obtain an up-to-date version.


Quick Installation and Test 
---------------------------

Follow these instruction to test i-PI. These assume to be run from a Linux 
environment, with a recent version of Python, Numpy and gfortran, and that 
the terminal is initially in the i-pi package directory (the directory 
containing this file).

* Generate the driver code

::

$ cd driver
$ make
$ cd ..

* Run one of the examples

This will first start the wrapper in the background, redirecting the output on 
a log file, then run a couple of instances of the driver code and then follow
the progress of the wrapper by monitoring the log file::

$ cd examples/tutorial/tutorial-1/
$ ../../../i-pi tutorial-1.xml > log &
$ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
$ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
$ tail -f log

The monitoring can be interrupted with ``CTRL+C`` when the run has finished (5000 steps)

