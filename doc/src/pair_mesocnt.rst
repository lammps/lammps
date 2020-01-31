.. index:: pair\_style mesocnt

pair\_style mesocnt command
===========================

Syntax
""""""


.. parsed-literal::

   pair_style mesocnt

Examples
""""""""


.. parsed-literal::

   pair_style mesocnt
   pair_coeff \* \* 10_10.cnt

Description
"""""""""""

Style *mesocnt* implements a mesoscopic potential
for the interaction of carbon nanotubes (CNTs). In this potential,
CNTs are modelled as chains of cylindrical segments in which
each infinitesimal surface element interacts with all other
CNT surface elements with the Lennard-Jones (LJ) term adopted from
the :doc:`airebo <pair_airebo>` style. The interaction energy
is then computed by integrating over the surfaces of all interacting
CNTs.

The potential is based on interactions between one cylindrical
segment and infinitely or semi-infinitely long CNTs as described
in :ref:`(Volkov1) <Volkov1>`. Chains of segments are
converted to these (semi-)infinite CNTs bases on an approximate
chain approach outlined in :ref:`(Volkov2) <Volkov2>`.
This allows to simplify the computation of the interactions
significantly and reduces the computational times to the
same order of magnitude as for regular bead spring models
where beads interact with the standard :doc:`pair_lj/cut <pair_lj>`
potential.

In LAMMPS, cylindrical segments are represented by bonds. Each
segment is defined by its two end points ("nodes") which correspond
to atoms in LAMMPS. For the exact functional form of the potential
and implementation details, the reader is referred to the 
original papers :ref:`(Volkov1) <Volkov1>` and 
:ref:`(Volkov2) <Volkov2>`.

The potential requires tabulated data provided in a single ASCII 
text file specified in the :doc:`pair_coeff <pair_coeff>` command. 
The first line of the file provides a time stamp and
general information. The second line lists four integers giving
the number of data points provided in the subsequent four
data tables. The third line lists four floating point numbers: 
the CNT radius R, the LJ parameter sigma and two numerical 
parameters delta1 and delta2. These four parameters are given
in Angstroms. This is followed by four data tables each separated
by a single empty line. The first two tables have two columns
and list the parameters uInfParallel and Gamma respectively.
The last two tables have three columns giving data on a quadratic
array and list the parameters Phi and uSemiParallel respectively.
uInfParallel and uSemiParallel are given in eV/Angstrom, Phi is
given in eV and Gamma is unitless.

Potential files for CNTs can be readily generated using the freely 
available code provided on

.. parsed-literal::
  
   https://github.com/phankl/cntpot

Using the same approach, it should also be possible to
generate potential files for other 1D systems such as
boron nitride nanotubes.

.. note::

   LAMMPS comes with one *mesocnt* style potential file
   where the default number of data points per table is 1001.
   This is sufficient for NVT simulations. For proper energy
   conservation, we recommend using a potential file where
   the resolution for Phi is at least 2001 data points.

.. note::

   The *mesocnt* style requires CNTs to be represented
   as a chain of atoms connected by bonds. Atoms need
   to be numbered consecutively within one chain. 
   Atoms belonging to different CNTs need to be assigned
   different molecule IDs.

A full summary of the method and LAMMPS implementation details
is expected to soon become available in Computer Physics
Communications.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support mixing.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

The *mesocnt* pair style do not write their information to :doc:`binary restart files <restart>`, 
since it is stored in tabulated potential files.
Thus, you need to re-specify the pair\_style and pair\_coeff commands in
an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This style is part of the USER-MISC package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair potential requires the :doc:`newton <newton>` setting to be
"on" for pair interactions.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------


.. _Volkov1:



**(Volkov1)** Volkov and Zhigilei, J Phys Chem C, 114, 5513 (2010).

.. _Volkov2:



**(Volkov2)** Volkov, Simov and Zhigilei, APS Meeting Abstracts, 
Q31.013 (2008).


