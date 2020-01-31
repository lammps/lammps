.. index:: temper/grem

temper/grem command
===================

Syntax
""""""


.. parsed-literal::

   temper/grem N M lambda fix-ID thermostat-ID seed1 seed2 index

* N = total # of timesteps to run
* M = attempt a tempering swap every this many steps
* lambda = initial lambda for this ensemble
* fix-ID = ID of *fix grem*
* thermostat-ID = ID of the thermostat that controls kinetic temperature
* seed1 = random # seed used to decide on adjacent temperature to partner with
* seed2 = random # seed for Boltzmann factor in Metropolis swap
* index = which temperature (0 to N-1) I am simulating (optional)

Examples
""""""""


.. parsed-literal::

   temper/grem 100000 1000 $\ *lambda* fxgREM fxnvt 0 58728
   temper/grem 40000 100 $\ *lambda* fxgREM fxnpt 0 32285 $\ *walkers*

Description
"""""""""""

Run a parallel tempering or replica exchange simulation in LAMMPS
partition mode using multiple generalized replicas (ensembles) of a
system defined by :doc:`fix grem <fix_grem>`, which stands for the
generalized replica exchange method (gREM) originally developed by
:ref:`(Kim) <KimStraub>`.  It uses non-Boltzmann ensembles to sample over first
order phase transitions. The is done by defining replicas with an
enthalpy dependent effective temperature

Two or more replicas must be used.  See the :doc:`temper <temper>`
command for an explanation of how to run replicas on multiple
partitions of one or more processors.

This command is a modification of the :doc:`temper <temper>` command and
has the same dependencies, restraints, and input variables which are
discussed there in greater detail.

Instead of temperature, this command performs replica exchanges in
lambda as per the generalized ensemble enforced by :doc:`fix grem <fix_grem>`.  The desired lambda is specified by *lambda*\ ,
which is typically a variable previously set in the input script, so
that each partition is assigned a different temperature.  See the
:doc:`variable <variable>` command for more details.  For example:


.. parsed-literal::

   variable lambda world 400 420 440 460
   fix fxnvt all nvt temp 300.0 300.0 100.0
   fix fxgREM all grem $\ *lambda* -0.05 -50000 fxnvt
   temper 100000 100 $\ *lambda* fxgREM fxnvt 3847 58382

would define 4 lambdas with constant kinetic temperature but unique
generalized temperature, and assign one of them to :doc:`fix grem <fix_grem>` used by each replica, and to the grem command.

As the gREM simulation runs for *N* timesteps, a swap between adjacent
ensembles will be attempted every *M* timesteps.  If *seed1* is 0,
then the swap attempts will alternate between odd and even pairings.
If *seed1* is non-zero then it is used as a seed in a random number
generator to randomly choose an odd or even pairing each time.  Each
attempted swap of temperatures is either accepted or rejected based on
a Metropolis criterion, derived for gREM by :ref:`(Kim) <Kim>`, which uses
*seed2* in the random number generator.

File management works identical to the :doc:`temper <temper>` command.
Dump files created by this fix contain continuous trajectories and
require post-processing to obtain per-replica information.

The last argument *index* in the grem command is optional and is used
when restarting a run from a set of restart files (one for each
replica) which had previously swapped to new lambda.  This is done
using a variable. For example if the log file listed the following for
a simulation with 5 replicas:


.. parsed-literal::

   500000 2 4 0 1 3

then a setting of


.. parsed-literal::

   variable walkers world 2 4 0 1 3

would be used to restart the run with a grem command like the example
above with $\ *walkers* as the last argument. This functionality is
identical to :doc:`temper <temper>`.


----------


Restrictions
""""""""""""


This command can only be used if LAMMPS was built with the USER-MISC
package.  See the :doc:`Build package <Build_package>` doc
page for more info.

This command must be used with :doc:`fix grem <fix_grem>`.

Related commands
""""""""""""""""

:doc:`fix grem <fix_grem>`, :doc:`temper <temper>`, :doc:`variable <variable>`

**Default:** none

.. _KimStraub:



**(Kim)** Kim, Keyes, Straub, J Chem Phys, 132, 224107 (2010).
