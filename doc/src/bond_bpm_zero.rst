.. index:: bond_style bpm/zero

bond_style bpm/zero command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/zero keyword value attribute1 attribute2 ...

* zero or more keyword/value pairs may be appended
* keyword = *manybody* or *overlay/pair* or *break* or *store/local* or *write/reference* or *read/reference*

  .. parsed-literal::

       *manybody* value = *yes* or *no*
          exercise the per-atom multibody communication machinery

       *overlay/pair* value = *yes* or *no*
          bonded particles will still interact with pair forces

       *break* value = *yes* or *no*
          indicates whether bonds break during a run

       *store/local* values = fix_ID N attributes ...
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

       *write/reference* values = fix_ID N
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

       *read/reference* value = filename
          (as for :doc:`bond_style bpm/spring <bond_bpm_spring>`)

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/zero
   bond_coeff 1 0.1

   bond_style bpm/zero break no
   bond_coeff 1 0.0

   bond_style bpm/zero write/reference myfix 1000 read/reference bond.ref
   dump 1 all local 1000 bond*.ref f_myfix[*]

Description
"""""""""""

.. versionadded:: TBD

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
*store/local*, *write/reference*, and *read/reference* keywords behave as
for :doc:`bond_style bpm/spring <bond_bpm_spring>`; see the :doc:`BPM Howto
<Howto_bpm>`.

----------

Restart and other info
"""""""""""""""""""""""""

This bond style writes the reference state of each bond to :doc:`binary
restart files <restart>`.  The reference state is not written to data
files.  After restarting, or after reading a data file, bonds reference
data can be restored using the *read/reference* option.  If *store/local*
is used, an internal fix records broken-bond data accessible through a
:doc:`dump local <dump>` command, as for the other BPM bond styles.  The
*write/reference* option likewise transfers the internal bond history to
an internal fix that can be output with a :doc:`dump local <dump>`
command; see the :doc:`BPM Howto <Howto_bpm>`.

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
