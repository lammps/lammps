.. index:: bond_style bpm/zero

bond_style bpm/zero command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/zero keyword value attribute1 attribute2 ...

* zero or more keyword/value pairs may be appended
* keyword = *manybody* or *overlay/pair* or *break* or *store/local* or *write/history* or *read/history*

  .. parsed-literal::

       *manybody* value = *yes* or *no*
          exercise the per-atom multibody communication machinery

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

       *store/local* values = fix_ID N attributes ...
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

       *write/history* values = fix_ID N
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

       *read/history* value = filename
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/zero
   bond_coeff 1 0.1

   bond_style bpm/zero break no
   bond_coeff 1 0.0

   bond_style bpm/zero write/history myfix 1000 read/history bond.ref
   dump 1 all local 1000 bond*.ref f_myfix[*]

Description
"""""""""""

.. versionadded:: 2Sep2026

The *bpm/zero* bond style is the :doc:`BPM package <Howto_bpm>` analogue of
:doc:`bond_style zero <bond_zero>`: it stores the initial reference state
of each bond and can break bonds individually, but it computes **no bond
force or energy**.  Like the other :doc:`BPM bond styles <bond_bpm_spring>`
the reference length is recorded when a bond is first computed in the setup
of a run, is preserved across run commands, and is written to :doc:`binary
restart files <restart>`.

A bond breaks when its strain :math:`(r - r_0)/r_0` exceeds the critical
value :math:`\epsilon_c` given by :doc:`bond_coeff <bond_coeff>` (unless
*break* is set to *no*).  Because no force is applied, broken or unbroken
the particles move only under the other forces in the system.

This style is intended for testing, debugging, and as a starting template,
not for production mechanics.

The following coefficient must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command:

* :math:`\epsilon_c`   (unitless), critical strain for breaking

The *manybody* keyword toggles an internal per-atom property and its
forward/reverse communication so the multibody machinery shared by the BPM
styles can be exercised in isolation.  The *overlay/pair*, *break*,
*store/local*, *write/history*, and *read/history* keywords behave as
for :doc:`bond_style bpm/spring <bond_bpm_spring>`; see the :doc:`BPM Howto
<Howto_bpm>`.

----------

Restart and other info
"""""""""""""""""""""""""


This bond style writes the history data (e.g. reference state) of each
bond and its per-type coefficients to :doc:`binary restart files <restart>`.
Loading a restart file restores bonds and their history data.  The history
data is NOT written to data files.  Reading a data file will therefore not
restore bond history such that bond reference states will be redefined.
Alternatively, bond history data can be saved and restored using the
*write/history* and *read/history* options.

If the *store/local* option is used, an internal fix will calculate
a local vector or local array depending on the number of input values.
The length of the vector or number of rows in the array is the number
of recorded, broken bonds.  If a single input is specified, a local
vector is produced. If two or more inputs are specified, a local array
is produced where the number of columns = the number of inputs.  The
vector or array can be accessed by any command that uses local values
from a compute as input. See the :doc:`Howto output <Howto_output>` page
for an overview of LAMMPS output options.

The vector or array will be floating point values that correspond to
the specified attribute.

Any settings with the *store/local* option are not saved to a restart
file and must be redefined.

If the *write/history* keyword is used, an internal fix will process
and transfer the internal bond history (e.g. :math:`r_0`) to an
internal fix labeled *fix_ID*.  This allows the internal bond history data,
as well as the IDs of the two atoms in the bond,  to be accessed by other
LAMMPS commands, in particular, :doc:`dump local <dump>`.

If the *read/history* keyword is used, history data is read from the file
labeled *filename*.  This is expected to be in the format of a LAMMPS
local dump file. The first two columns of the history file must contain
the IDs of the two atoms in the bond, with the remaining columns corresponding
 to the internal bond data.  For more details on formatting the remaining file
 see :doc:`Howto bpm <Howto_bpm>`.

Restrictions
""""""""""""

This bond style is part of the BPM package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

As with the other BPM bond styles, :doc:`newton <newton>` must be set to
*bond off* and the special bond weights must be

.. code-block:: LAMMPS

   special_bonds lj 0 1 1 coul 1 1 1

(or all weights one with *overlay/pair yes*).

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`bond_style zero <bond_zero>`,
:doc:`bond_style bpm/spring <bond_bpm_spring>`, :doc:`Howto BPM <Howto_bpm>`

Default
"""""""

The option default is *manybody* = *no* (plus the BPM defaults
*overlay/pair* = *no*, *break* = *yes*).
