.. index:: pair_style tracker

pair_style tracker command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style tracker keyword

* zero or more keyword/arg pairs may be appended
* keyword = *finite*

  .. parsed-literal::

      *finite* value = none
         pair style uses atomic diameters to identify contacts

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay tracker ...
   pair_coeff 1 1 tracker 2.0

   pair_style hybrid/overlay tracker finite ...
   pair_coeff * * tracker

   fix 1 all pair/tracker 1000 time/created time/broken
   dump 1 all local 1000 dump.local f_1[1] f_1[2]
   dump_modify 1 write_header no

Description
"""""""""""

Style *tracker* monitors information about pairwise interactions.
It does not calculate any forces on atoms.
:doc:`Pair hybrid/overlay <pair_hybrid>` can be used to combine this pair
style with another pair style. Style *tracker*  must be used in conjunction
with about :doc:`fix pair_tracker <fix_pair_tracker>` which contains
information on what data can be output.

If the *finite* keyword is not defined, the following coefficients must be
defined for each pair of atom types via the :doc:`pair_coeff <pair_coeff>`
command as in the examples above, or in the data file or restart files
read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutoff (distance units)

If the *finite* keyword is defined, no coefficients may be defined.
Interaction cutoffs are alternatively calculated based on the
diameter of finite particles.


Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the cutoff coefficient and cutoff
distance for this pair style can be mixed.  The cutoff is always mixed via a
*geometric* rule.  The cutoff is mixed according to the pair_modify
mix value.  The default mix value is *geometric*\ .  See the
"pair_modify" command for details.

This pair style writes its information to :doc:`binary restart files <restart>`, so
pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

The :doc:`pair_modify <pair_modify>` shift, table, and tail options
are not relevant for this pair style.

----------

Restrictions
""""""""""""

A corresponding :doc:`fix pair_tracker <fix_pair_tracker>` must be defined
to use this pair style.

This pair style is currently incompatible with granular pair styles that extend
beyond the contact (e.g. JKR and DMT).

This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix pair_tracker <fix_pair_tracker>`

Default
"""""""

none
