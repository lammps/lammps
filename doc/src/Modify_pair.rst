Pair styles
===========

Classes that compute pairwise interactions are derived from the Pair
class.  In LAMMPS, pairwise calculation include many-body potentials
such as EAM or Tersoff where particles interact without a static bond
topology.  New styles can be created to add new pair potentials to
LAMMPS.

Pair\_lj\_cut.cpp is a simple example of a Pair class, though it
includes some optional methods to enable its use with rRESPA.

Here is a brief description of the class methods in pair.h:

+---------------------------------+-------------------------------------------------------------------+
| compute                         | workhorse routine that computes pairwise interactions             |
+---------------------------------+-------------------------------------------------------------------+
| settings                        | reads the input script line with arguments you define             |
+---------------------------------+-------------------------------------------------------------------+
| coeff                           | set coefficients for one i,j type pair                            |
+---------------------------------+-------------------------------------------------------------------+
| init\_one                       | perform initialization for one i,j type pair                      |
+---------------------------------+-------------------------------------------------------------------+
| init\_style                     | initialization specific to this pair style                        |
+---------------------------------+-------------------------------------------------------------------+
| write & read\_restart           | write/read i,j pair coeffs to restart files                       |
+---------------------------------+-------------------------------------------------------------------+
| write & read\_restart\_settings | write/read global settings to restart files                       |
+---------------------------------+-------------------------------------------------------------------+
| single                          | force and energy of a single pairwise interaction between 2 atoms |
+---------------------------------+-------------------------------------------------------------------+
| compute\_inner/middle/outer     | versions of compute used by rRESPA                                |
+---------------------------------+-------------------------------------------------------------------+

The inner/middle/outer routines are optional.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
