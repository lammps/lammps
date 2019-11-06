.. index:: fix mscg

fix mscg command
================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID mscg N keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mscg = style name of this fix command
* N = envoke this fix every this many timesteps
* zero or more keyword/value pairs may be appended
* keyword = *range* or *name* or *max*
  
  .. parsed-literal::
  
       *range* arg = *on* or *off*
         *on* = range finding functionality is performed
         *off* = force matching functionality is performed
       *name* args = name1 ... nameN
         name1,...,nameN = string names for each atom type (1-Ntype)
       *max* args = maxb maxa maxd
         maxb,maxa,maxd = maximum bonds/angles/dihedrals per atom



Examples
""""""""


.. parsed-literal::

   fix 1 all mscg 1
   fix 1 all mscg 1 range name A B
   fix 1 all mscg 1 max 4 8 20

Description
"""""""""""

This fix applies the Multi-Scale Coarse-Graining (MSCG) method to
snapshots from a dump file to generate potentials for coarse-grained
simulations from all-atom simulations, using a force-matching
technique (:ref:`Izvekov <Izvekov>`, :ref:`Noid <Noid>`).

It makes use of the MS-CG library, written and maintained by Greg
Voth's group at the University of Chicago, which is freely available
on their `MS-CG GitHub site <https://github.com/uchicago-voth/MSCG-release>`_.  See instructions
on obtaining and installing the MS-CG library in the src/MSCG/README
file, which must be done before you build LAMMPS with this fix command
and use the command in a LAMMPS input script.

An example script using this fix is provided the examples/mscg
directory.

The general workflow for using LAMMPS in conjunction with the MS-CG
library to create a coarse-grained model and run coarse-grained
simulations is as follows:

1. Perform all-atom simulations on the system to be coarse grained.
2. Generate a trajectory mapped to the coarse-grained model.
3. Create input files for the MS-CG library.
4. Run the range finder functionality of the MS-CG library.
5. Run the force matching functionality of the MS-CG library.
6. Check the results of the force matching.
7. Run coarse-grained simulations using the new coarse-grained potentials.

This fix can perform the range finding and force matching steps 4 and
5 of the above workflow when used in conjunction with the
:doc:`rerun <rerun>` command.  It does not perform steps 1-3 and 6-7.

Step 2 can be performed using a Python script (what is the name?)
provided with the MS-CG library which defines the coarse-grained model
and converts a standard LAMMPS dump file for an all-atom simulation
(step 1) into a LAMMPS dump file which has the positions of and forces
on the coarse-grained beads.

In step 3, an input file named "control.in" is needed by the MS-CG
library which sets parameters for the range finding and force matching
functionalities.  See the examples/mscg/control.in file as an example.
And see the documentation provided with the MS-CG library for more
info on this file.

When this fix is used to perform steps 4 and 5, the MS-CG library also
produces additional output files.  The range finder functionality
(step 4) outputs files defining pair and bonded interaction ranges.
The force matching functionality (step 5) outputs tabulated force
files for every interaction in the system. Other diagnostic files can
also be output depending on the parameters in the MS-CG library input
script.  Again, see the documentation provided with the MS-CG library
for more info.


----------


The *range* keyword specifies which MS-CG library functionality should
be invoked. If *on*\ , the step 4 range finder functionality is invoked.
*off*\ , the step 5 force matching functionality is invoked.

If the *name* keyword is used, string names are defined to associate
with the integer atom types in LAMMPS.  *Ntype* names must be
provided, one for each atom type (1-Ntype).

The *max* keyword specifies the maximum number of bonds, angles, and
dihedrals a bead can have in the coarse-grained model.

Restrictions
""""""""""""


This fix is part of the MSCG package. It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

The MS-CG library uses C++11, which may not be supported by older
compilers. The MS-CG library also has some additional numeric library
dependencies, which are described in its documentation.

Currently, the MS-CG library is not setup to run in parallel with MPI,
so this fix can only be used in a serial LAMMPS build and run
on a single processor.

**Related commands:** none

Default
"""""""

The default keyword settings are range off, max 4 12 36.


----------


.. _Izvekov:



**(Izvekov)** Izvekov, Voth, J Chem Phys 123, 134105 (2005).

.. _Noid:



**(Noid)** Noid, Chu, Ayton, Krishna, Izvekov, Voth, Das, Andersen, J
Chem Phys 128, 134105 (2008).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
