.. index:: fix reax/c/bonds

fix reax/c/bonds command
========================

fix reax/c/bonds/kk command
===========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID reaxc/bonds Nevery filename

* ID, group-ID are documented in :doc:`fix <fix>` command
* reax/bonds = style name of this fix command
* Nevery = output interval in timesteps
* filename = name of output file

Examples
""""""""


.. parsed-literal::

   fix 1 all reax/c/bonds 100 bonds.reaxc

Description
"""""""""""

Write out the bond information computed by the ReaxFF potential specified
by :doc:`pair_style reax/c <pair_reaxc>` in the exact same format as the
original stand-alone ReaxFF code of Adri van Duin.  The bond information
is written to *filename* on timesteps that are multiples of *Nevery*\ ,
including timestep 0.  For time-averaged chemical species analysis,
please see the :doc:`fix reaxc/c/species <fix_reaxc_species>` command.

The specified group-ID is ignored by this fix.

The format of the output file should be reasonably self-explanatory.
The meaning of the column header abbreviations is as follows:

* id = atom id
* type = atom type
* nb = number of bonds
* id\_1 = atom id of first bond
* id\_nb = atom id of Nth bond
* mol = molecule id
* bo\_1 = bond order of first bond
* bo\_nb = bond order of Nth bond
* abo = atom bond order (sum of all bonds)
* nlp = number of lone pairs
* q = atomic charge

If the filename ends with ".gz", the output file is written in gzipped
format.  A gzipped dump file will be about 3x smaller than the text
version, but will also take longer to write.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed in :doc:`Speed <Speed>`
of the manual.  The accelerated styles take the same arguments and
should produce the same results, except for round-off and precision
issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See :doc:`Speed <Speed>` of the manual for
more instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


The fix reax/c/bonds command requires that the :doc:`pair_style reax/c <pair_reaxc>` is invoked.  This fix is part of the
USER-REAXC package.  It is only enabled if LAMMPS was built with that
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

To write gzipped bond files, you must compile LAMMPS with the
-DLAMMPS\_GZIP option.

Related commands
""""""""""""""""

:doc:`pair_style reax/c <pair_reaxc>`, :doc:`fix reax/c/species <fix_reaxc_species>`

**Default:** none
