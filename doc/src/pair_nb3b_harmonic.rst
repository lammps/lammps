.. index:: pair_style nb3b/harmonic

pair_style nb3b/harmonic command
================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style nb3b/harmonic

Examples
""""""""

.. code-block:: LAMMPS

   pair_style nb3b/harmonic
   pair_coeff * * MgOH.nb3bharmonic Mg O H

Description
"""""""""""

This pair style computes a non-bonded 3-body harmonic potential for the
energy E of a system of atoms as

.. math::

   E = K (\theta - \theta_0)^2

where :math:`\theta_0` is the equilibrium value of the angle and *K* is a
prefactor. Note that the usual 1/2 factor is included in *K*\ . The form
of the potential is identical to that used in angle_style *harmonic*\ ,
but in this case, the atoms do not need to be explicitly bonded.

Only a single pair_coeff command is used with this style which
specifies a potential file with parameters for specified elements.
These are mapped to LAMMPS atom types by specifying N additional
arguments after the filename in the pair_coeff command, where N is the
number of LAMMPS atom types:

* filename
* N element names = mapping of elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine a file SiC.nb3b.harmonic has potential values
for Si and C.  If your LAMMPS simulation has 4 atoms types and you
want the first 3 to be Si, and the fourth to be C, you would use the
following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * SiC.nb3b.harmonic Si Si Si C

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the potential file.  The final C argument maps LAMMPS atom
type 4 to the C element in the potential file.  If a mapping value is
specified as NULL, the mapping is not performed.  This can be used
when the potential is used as part of the *hybrid* pair style.  The
NULL values are placeholders for atom types that will be used with
other potentials. An example of a pair_coeff command for use with the
*hybrid* pair style is:

.. code-block:: LAMMPS

   pair_coeff * * nb3b/harmonic MgOH.nb3b.harmonic Mg O H

Three-body non-bonded harmonic files in the *potentials* directory of
the LAMMPS distribution have a ".nb3b.harmonic" suffix.  Lines that
are not blank or comments (starting with #) define parameters for a
triplet of elements.

Each entry has six arguments. The first three are atom types as
referenced in the LAMMPS input file. The first argument specifies the
central atom. The fourth argument indicates the *K* parameter. The
fifth argument indicates :math:`\theta_0`. The sixth argument indicates a
separation cutoff in Angstroms.

For a given entry, if the second and third arguments are identical,
then the entry is for a cutoff for the distance between types 1 and 2
(values for *K* and :math:`\theta_0` are irrelevant in this case).

For a given entry, if the first three arguments are all different, then
the entry is for the *K* and :math:`\theta_0` parameters (the cutoff in
this case is irrelevant).

It is required that the potential file contains entries for *all*
permutations of the elements listed in the pair_coeff command.
If certain combinations are not parameterized the corresponding
parameters should be set to zero. The potential file can also
contain entries for additional elements which are not used in
a particular simulation; LAMMPS ignores those entries.

----------

Restrictions
""""""""""""

This pair style can only be used if LAMMPS was built with the MANYBODY
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none
