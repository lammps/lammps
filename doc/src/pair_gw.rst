.. index:: pair_style gw

pair_style gw command
=====================

pair_style gw/zbl command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *gw* or *gw/zbl*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style gw
   pair_coeff * * SiC.gw Si C C

   pair_style gw/zbl
   pair_coeff * * SiC.gw.zbl C Si

Description
"""""""""""

The *gw* style computes a 3-body :ref:`Gao-Weber <Gao>` potential;
similarly *gw/zbl* combines this potential with a modified
repulsive ZBL core function in a similar fashion as implemented
in the :doc:`tersoff/zbl <pair_tersoff_zbl>` pair style.

Unfortunately the author of this contributed code has not been
able to submit a suitable documentation explaining the details
of the potentials. The LAMMPS developers thus have finally decided
to release the code anyway with only the technical explanations.
For details of the model and the parameters, please refer to the
linked publication.

Only a single pair_coeff command is used with the *gw* and *gw/zbl*
styles which specifies a Gao-Weber potential file with parameters
for all needed elements.  These are mapped to LAMMPS atom types by
specifying N additional arguments after the filename in the pair_coeff
command, where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of GW elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine a file SiC.gw has Gao-Weber values for Si and C.
If your LAMMPS simulation has 4 atoms types and you want the first 3 to
be Si, and the fourth to be C, you would use the following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * SiC.gw Si Si Si C

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the GW file.  The final C argument maps LAMMPS atom type 4
to the C element in the GW file.  If a mapping value is specified as
NULL, the mapping is not performed.  This can be used when a *gw*
potential is used as part of the *hybrid* pair style.  The NULL values
are placeholders for atom types that will be used with other
potentials.

Gao-Weber files in the *potentials* directory of the LAMMPS
distribution have a ".gw" suffix.  Gao-Weber with ZBL files
have a ".gz.zbl" suffix. The structure of the potential files
is similar to other many-body potentials supported by LAMMPS.
You have to refer to the comments in the files and the literature
to learn more details.

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above from values in the potential file.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package. It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Gao-Weber potential files provided with LAMMPS (see the
potentials directory) are parameterized for metal :doc:`units <units>`.
You can use the GW potential with any LAMMPS units, but you would need
to create your own GW potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none

----------

.. _Gao:

**(Gao)** Gao and Weber, Nuclear Instruments and Methods in Physics
Research B 191 (2012) 504.
