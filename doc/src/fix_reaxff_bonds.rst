.. index:: fix reaxff/bonds
.. index:: fix reaxff/bonds/kk

fix reaxff/bonds command
========================

Accelerator Variants: *reaxff/bonds/kk*

Syntax
""""""

.. parsed-literal::

   fix ID group-ID reaxff/bonds Nevery filename

* ID, group-ID are documented in :doc:`fix <fix>` command
* reax/bonds = style name of this fix command
* Nevery = output interval in timesteps
* filename = name of output file

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all reaxff/bonds 100 bonds.reaxff

Description
"""""""""""

Write out the bond information computed by the ReaxFF potential specified
by :doc:`pair_style reaxff <pair_reaxff>` in the exact same format as the
original stand-alone ReaxFF code of Adri van Duin.  The bond information
is written to *filename* on timesteps that are multiples of *Nevery*,
including timestep 0.  For time-averaged chemical species analysis,
please see the :doc:`fix reaxff/species <fix_reaxff_species>` command.

The specified group-ID is ignored by this fix.

The format of the output file should be reasonably self-explanatory.
The meaning of the column header abbreviations is as follows:

* id = atom id
* type = atom type
* nb = number of bonds
* id_1 = atom id of first bond
* id_nb = atom id of Nth bond
* mol = molecule id
* bo_1 = bond order of first bond
* bo_nb = bond order of Nth bond
* abo = atom bond order (sum of all bonds)
* nlp = number of lone pairs
* q = atomic charge

If the filename ends with ".gz", the output file is written in gzipped
format.  A gzipped dump file will be about 3x smaller than the text
version, but will also take longer to write.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The fix reaxff/bonds command requires that the :doc:`pair_style reaxff <pair_reaxff>` is invoked.  This fix is part of the
REAXFF package.  It is only enabled if LAMMPS was built with that
package.  See the :doc:`Build package <Build_package>` page for more
info.

To write gzipped bond files, you must compile LAMMPS with the
-DLAMMPS_GZIP option.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix reaxff/species <fix_reaxff_species>`

Default
"""""""

none
