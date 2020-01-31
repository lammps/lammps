.. index:: pair\_style lcbop

pair\_style lcbop command
=========================

Syntax
""""""


.. parsed-literal::

   pair_style lcbop

Examples
""""""""


.. parsed-literal::

   pair_style lcbop
   pair_coeff \* \* ../potentials/C.lcbop C

Description
"""""""""""

The *lcbop* pair style computes the long-range bond-order potential
for carbon (LCBOP) of :ref:`(Los and Fasolino) <Los>`.  See section II in
that paper for the analytic equations associated with the potential.

Only a single pair\_coeff command is used with the *lcbop* style which
specifies an LCBOP potential file with parameters for specific
elements.  These are mapped to LAMMPS atom types by specifying N
additional arguments after the filename in the pair\_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of LCBOP elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, if your LAMMPS simulation has 4 atom types and you want
the 1st 3 to be C you would use the following pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* C.lcbop C C C NULL

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first C argument maps LAMMPS atom type 1 to the C element in the
LCBOP file. If a mapping value is specified as NULL, the mapping is
not performed.  This can be used when a *lcbop* potential is used as
part of the *hybrid* pair style.  The NULL values are placeholders for
atom types that will be used with other potentials.

The parameters/coefficients for the LCBOP potential as applied to C
are listed in the C.lcbop file to agree with the original :ref:`(Los and Fasolino) <Los>` paper.  Thus the parameters are specific to this
potential and the way it was fit, so modifying the file should be done
carefully.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

Restrictions
""""""""""""


This pair styles is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair potential requires the :doc:`newton <newton>` setting to be
"on" for pair interactions.

The C.lcbop potential file provided with LAMMPS (see the potentials
directory) is parameterized for metal :doc:`units <units>`.  You can use
the LCBOP potential with any LAMMPS units, but you would need to
create your own LCBOP potential file with coefficients listed in the
appropriate units if your simulation doesn't use "metal" units.

Related commands
""""""""""""""""

:doc:`pair_airebo <pair_airebo>`, :doc:`pair_coeff <pair_coeff>`

**Default:** none


----------


.. _Los:



**(Los and Fasolino)** J. H. Los and A. Fasolino, Phys. Rev. B 68, 024107
(2003).


