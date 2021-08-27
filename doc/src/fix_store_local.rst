.. index:: fix store/local

fix store/local command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID store/local N nvalues

* ID, group-ID are documented in :doc:`fix <fix>` command
* store/local = style name of this fix command
* N = prepare data for output every this many timesteps
* nvalues = number of values stored by this fix

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all store/local 1000 2
   dump 1 all local 1000 dump.local f_1[1] f_1[2]

Description
"""""""""""

This fix provides the ability to store local data produced by
some LAMMPS commands including some pair and bond styles so it can be output.
Data is accumulated over a span of *N* timesteps before being deleted.
The number of datums generated, aggregated across all processors, depends on
the associated commands. Data is only included if it is generated from atoms
within the fix group-ID.

----------

Restart, fix_modify, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Accumulated local data is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Output info
"""""""""""

This compute calculates a local vector or local array depending on the
number of input values.  The length of the vector or number of rows in
the array is the number of recorded, lost interactions.  If a single input is
specified, a local vector is produced.  If two or more inputs are
specified, a local array is produced where the number of columns = the
number of inputs.  The vector or array can be accessed by any command
that uses local values from a compute as input.  See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The vector or array values will be doubles that correspond to the
specified attribute.

Restrictions
""""""""""""

Must be used in conjunction with another LAMMPS class which outputs local data.

Related commands
""""""""""""""""

:doc:`pair tracker <pair_tracker>`
:doc:`bond bpm/rotational <bond_bpm_rotational>`
:doc:`dump local <dump>`

Default
"""""""

none
