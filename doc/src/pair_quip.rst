.. index:: pair_style quip

pair_style quip command
=======================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style quip

Examples
""""""""


.. code-block:: LAMMPS

   pair_style      quip
   pair_coeff      * * gap_example.xml "Potential xml_label=GAP_2014_5_8_60_17_10_38_466" 14
   pair_coeff      * * sw_example.xml "IP SW" 14

Description
"""""""""""

Style *quip* provides an interface for calling potential routines from
the QUIP package. QUIP is built separately, and then linked to
LAMMPS. The most recent version of the QUIP package can be downloaded
from GitHub:
`https://github.com/libAtoms/QUIP <https://github.com/libAtoms/QUIP>`_. The
interface is chiefly intended to be used to run Gaussian Approximation
Potentials (GAP), which are described in the following publications:
:ref:`(Bartok et al) <Bartok_2010>` and :ref:`(PhD thesis of Bartok) <Bartok_PhD>`.

Only a single pair\_coeff command is used with the *quip* style that
specifies a QUIP potential file containing the parameters of the
potential for all needed elements in XML format. This is followed by a
QUIP initialization string. Finally, the QUIP elements are mapped to
LAMMPS atom types by specifying N atomic numbers, where N is the
number of LAMMPS atom types:

* QUIP filename
* QUIP initialization string
* N atomic numbers = mapping of QUIP elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

A QUIP potential is fully specified by the filename which contains the
parameters of the potential in XML format, the initialization string,
and the map of atomic numbers.

GAP potentials can be obtained from the Data repository section of
`http://www.libatoms.org <http://www.libatoms.org>`_, where the
appropriate initialization strings are also advised. The list of
atomic numbers must be matched to the LAMMPS atom types specified in
the LAMMPS data file or elsewhere.

Two examples input scripts are provided in the examples/USER/quip
directory.

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


This pair style is part of the USER-QUIP package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

QUIP potentials are parameterized in electron-volts and Angstroms and
therefore should be used with LAMMPS metal :doc:`units <units>`.

QUIP potentials are generally not designed to work with the scaling
factors set by the :doc:`special_bonds <special_bonds>` command.  The
recommended setting in molecular systems is to include all
interactions, i.e. to use *special\_bonds lj/coul 1.0 1.0 1.0*. Scaling
factors > 0.0 will be ignored and treated as 1.0. The only exception
to this rule is if you know that your QUIP potential needs to exclude
bonded, 1-3, or 1-4 interactions and does not already do this exclusion
within QUIP. Then a factor 0.0 needs to be used which will remove such
pairs from the neighbor list. This needs to be very carefully tested,
because it may remove pairs from the neighbor list that are still
required.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`


----------


.. _Bartok\_2010:



**(Bartok\_2010)** AP Bartok, MC Payne, R Kondor, and G Csanyi, Physical
Review Letters 104, 136403 (2010).

.. _Bartok\_PhD:



**(Bartok\_PhD)** A Bartok-Partay, PhD Thesis, University of Cambridge,
(2010).
