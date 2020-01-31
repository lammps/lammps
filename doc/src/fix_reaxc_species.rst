.. index:: fix reax/c/species

fix reax/c/species command
==========================

fix reax/c/species/kk command
=============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID reax/c/species Nevery Nrepeat Nfreq filename keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* reax/c/species = style name of this command
* Nevery = sample bond-order every this many timesteps
* Nrepeat = # of bond-order samples used for calculating averages
* Nfreq = calculate average bond-order every this many timesteps
* filename = name of output file
* zero or more keyword/value pairs may be appended
* keyword = *cutoff* or *element* or *position*
  
  .. parsed-literal::
  
       *cutoff* value = I J Cutoff
         I, J = atom types
         Cutoff = Bond-order cutoff value for this pair of atom types
       *element* value = Element1, Element2, ...
       *position* value = posfreq filepos
         posfreq = write position files every this many timestep
         filepos = name of position output file



Examples
""""""""


.. parsed-literal::

   fix 1 all reax/c/species 10 10 100 species.out
   fix 1 all reax/c/species 1 2 20 species.out cutoff 1 1 0.40 cutoff 1 2 0.55
   fix 1 all reax/c/species 1 100 100 species.out element Au O H position 1000 AuOH.pos

Description
"""""""""""

Write out the chemical species information computed by the ReaxFF
potential specified by :doc:`pair_style reax/c <pair_reaxc>`.
Bond-order values (either averaged or instantaneous, depending on
value of *Nrepeat*\ ) are used to determine chemical bonds.  Every
*Nfreq* timesteps, chemical species information is written to
*filename* as a two line output.  The first line is a header
containing labels. The second line consists of the following:
timestep, total number of molecules, total number of distinct species,
number of molecules of each species.  In this context, "species" means
a unique molecule.  The chemical formula of each species is given in
the first line.

If the filename ends with ".gz", the output file is written in gzipped
format.  A gzipped dump file will be about 3x smaller than the text version,
but will also take longer to write.

Optional keyword *cutoff* can be assigned to change the minimum
bond-order values used in identifying chemical bonds between pairs of
atoms.  Bond-order cutoffs should be carefully chosen, as bond-order
cutoffs that are too small may include too many bonds (which will
result in an error), while cutoffs that are too large will result in
fragmented molecules.  The default cutoff of 0.3 usually gives good
results.

The optional keyword *element* can be used to specify the chemical
symbol printed for each LAMMPS atom type. The number of symbols must
match the number of LAMMPS atom types and each symbol must consist of
1 or 2 alphanumeric characters. Normally, these symbols should be
chosen to match the chemical identity of each LAMMPS atom type, as
specified using the :doc:`reax/c pair\_coeff <pair_reaxc>` command and
the ReaxFF force field file.

The optional keyword *position* writes center-of-mass positions of
each identified molecules to file *filepos* every *posfreq* timesteps.
The first line contains information on timestep, total number of
molecules, total number of distinct species, and box dimensions.  The
second line is a header containing labels.  From the third line
downward, each molecule writes a line of output containing the
following information: molecule ID, number of atoms in this molecule,
chemical formula, total charge, and center-of-mass xyz positions of
this molecule.  The xyz positions are in fractional coordinates
relative to the box dimensions.

For the keyword *position*\ , the *filepos* is the name of the output
file.  It can contain the wildcard character "\*".  If the "\*"
character appears in *filepos*\ , then one file per snapshot is written
at *posfreq* and the "\*" character is replaced with the timestep
value.  For example, AuO.pos.\* becomes AuO.pos.0, AuO.pos.1000, etc.


----------


The *Nevery*\ , *Nrepeat*\ , and *Nfreq* arguments specify on what
timesteps the bond-order values are sampled to get the average bond
order.  The species analysis is performed using the average bond-order
on timesteps that are a multiple of *Nfreq*\ .  The average is over
*Nrepeat* bond-order samples, computed in the preceding portion of the
simulation every *Nevery* timesteps.  *Nfreq* must be a multiple of
*Nevery* and *Nevery* must be non-zero even if *Nrepeat* is 1.
Also, the timesteps
contributing to the average bond-order cannot overlap,
i.e. Nrepeat\*Nevery can not exceed Nfreq.

For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then bond-order
values on timesteps 90,92,94,96,98,100 will be used to compute the
average bond-order for the species analysis output on timestep 100.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes both a global vector of length 2 and a per-atom
vector, either of which can be accessed by various :doc:`output commands <Howto_output>`.  The values in the global vector are
"intensive".

The 2 values in the global vector are as follows:

* 1 = total number of molecules
* 2 = total number of distinct species

The per-atom vector stores the molecule ID for each atom as identified
by the fix.  If an atom is not in a molecule, its ID will be 0.
For atoms in the same molecule, the molecule ID for all of them
will be the same and will be equal to the smallest atom ID of
any atom in the molecule.

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


The "fix reax/c/species" currently only works with :doc:`pair_style reax/c <pair_reaxc>` and it requires that the :doc:`pair_style reax/c <pair_reaxc>` be invoked.  This fix is part of the
USER-REAXC package.  It is only enabled if LAMMPS was built with that
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

To write gzipped species files, you must compile LAMMPS with the
-DLAMMPS\_GZIP option.

It should be possible to extend it to other reactive pair\_styles (such as
:doc:`rebo <pair_airebo>`, :doc:`airebo <pair_airebo>`,
:doc:`comb <pair_comb>`, and :doc:`bop <pair_bop>`), but this has not yet been done.

Related commands
""""""""""""""""

:doc:`pair_style reax/c <pair_reaxc>`, :doc:`fix reax/c/bonds <fix_reaxc_bonds>`

Default
"""""""

The default values for bond-order cutoffs are 0.3 for all I-J pairs.  The
default element symbols are C, H, O, N.  Position files are not written
by default.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
