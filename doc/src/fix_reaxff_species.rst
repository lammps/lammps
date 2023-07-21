.. index:: fix reaxff/species
.. index:: fix reaxff/species/kk

fix reaxff/species command
==========================

Accelerator Variants: *reaxff/species/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID reaxff/species Nevery Nrepeat Nfreq filename keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* reaxff/species = style name of this command
* Nevery = sample bond-order every this many timesteps
* Nrepeat = # of bond-order samples used for calculating averages
* Nfreq = calculate average bond-order every this many timesteps
* filename = name of output file
* zero or more keyword/value pairs may be appended
* keyword = *cutoff* or *element* or *position* or *delete*

  .. parsed-literal::

       *cutoff* value = I J Cutoff
         I, J = atom types (see asterisk form below)
         Cutoff = Bond-order cutoff value for this pair of atom types
       *element* value = Element1, Element2, ...
       *position* value = posfreq filepos
         posfreq = write position files every this many timestep
         filepos = name of position output file
       *delete* value = filedel keyword value
         filedel = name of delete species output file
         keyword = *specieslist* or *masslimit*
           *specieslist* value = Nspecies Species1 Species2 ...
             Nspecies = number of species in list
           *masslimit* value = massmin massmax
             massmin = minimum molecular weight of species to delete
             massmax = maximum molecular weight of species to delete
       *delete_rate_limit* value = Nlimit Nsteps
             Nlimit = maximum number of deletions allowed to occur within interval
             Nsteps = the interval (number of timesteps) over which to count deletions

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all reaxff/species 10 10 100 species.out
   fix 1 all reaxff/species 1 2 20 species.out cutoff 1 1 0.40 cutoff 1 2*3 0.55
   fix 1 all reaxff/species 1 100 100 species.out element Au O H position 1000 AuOH.pos
   fix 1 all reaxff/species 1 100 100 species.out delete species.del masslimit 0 50

Description
"""""""""""

Write out the chemical species information computed by the ReaxFF
potential specified by :doc:`pair_style reaxff <pair_reaxff>`.
Bond-order values (either averaged or instantaneous, depending on
value of *Nrepeat*\ ) are used to determine chemical bonds.  Every
*Nfreq* timesteps, chemical species information is written to
*filename* as a two line output.  The first line is a header
containing labels. The second line consists of the following:
timestep, total number of molecules, total number of distinct species,
number of molecules of each species.  In this context, "species" means
a unique molecule.  The chemical formula of each species is given in
the first line.

.. warning::

   In order to compute averaged data, it is required that there are no
   neighbor list rebuilds for at least Nrepeat\*Nevery steps preceding
   each *Nfreq* step.  For that reason, fix *reaxff/species* may
   change your neighbor list settings.  Reneighboring will occur no
   more frequently than every Nrepeat\*Nevery timesteps, and will
   occur less frequently if *Nfreq* is not a multiple of
   Nrepeat\*Nevery.  There will be a warning message showing the new
   settings.  Having a *Nfreq* setting that is larger than what is
   required for correct computation of the ReaxFF force field
   interactions, in combination with certain *Nrepeat* and *Nevery*
   settings, can thus lead to incorrect results.  For typical ReaxFF
   calculations, reneighboring only every 100 steps is already quite a
   low frequency.

If the filename ends with ".gz", the output file is written in gzipped
format.  A gzipped dump file will be about 3x smaller than the text version,
but will also take longer to write.

.. versionadded:: 15Jun2023

   Support for wildcards added

Optional keyword *cutoff* can be assigned to change the minimum
bond-order values used in identifying chemical bonds between pairs of
atoms.  Bond-order cutoffs should be carefully chosen, as bond-order
cutoffs that are too small may include too many bonds (which will result
in an error), while cutoffs that are too large will result in fragmented
molecules.  The default cutoff of 0.3 usually gives good results.  A
wildcard asterisk can be used in place of or in conjunction with the I,J
arguments to set the bond-order cutoff for multiple pairs of atom types.
This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If :math:`N` is
the number of atom types, then an asterisk with no numeric values means
all types from 1 to :math:`N`.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from n to
:math:`N` (inclusive).  A middle asterisk means all types from m to n
(inclusive).

The optional keyword *element* can be used to specify the chemical
symbol printed for each LAMMPS atom type. The number of symbols must
match the number of LAMMPS atom types and each symbol must consist of
1 or 2 alphanumeric characters. Normally, these symbols should be
chosen to match the chemical identity of each LAMMPS atom type, as
specified using the :doc:`reaxff pair_coeff <pair_reaxff>` command and
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

For the keyword *position*, the *filepos* is the name of the output
file.  It can contain the wildcard character "\*".  If the "\*"
character appears in *filepos*, then one file per snapshot is written
at *posfreq* and the "\*" character is replaced with the timestep
value.  For example, AuO.pos.\* becomes AuO.pos.0, AuO.pos.1000, etc.

.. versionadded:: 3Aug2022

The optional keyword *delete* enables the periodic removal of
molecules from the system.  Criteria for deletion can be either a list
of specific chemical formulae or a range of molecular weights.
Molecules are deleted every *Nfreq* timesteps, and bond connectivity
is determined using the *Nevery* and *Nrepeat* keywords.  The
*filedel* argument is the name of the output file that records the
species that are removed from the system.  The *specieslist* keyword
permits specific chemical species to be deleted.  The *Nspecies*
argument specifies how many species are eligible for deletion and is
followed by a list of chemical formulae, whose strings are compared to
species identified by this fix.  For example, "specieslist 2 CO CO2"
deletes molecules that are identified as "CO" and "CO2" in the species
output file.  When using the *specieslist* keyword, the *filedel* file
has the following format: the first line lists the chemical formulae
eligible for deletion, and each additional line contains the timestep
on which a molecule deletion occurs and the number of each species
deleted on that timestep.  The *masslimit* keyword permits deletion of
molecules with molecular weights between *massmin* and *massmax*.
When using the *masslimit* keyword, each line of the *filedel* file
contains the timestep on which deletions occurs, followed by how many
of each species are deleted (with quantities preceding chemical
formulae).  The *specieslist* and *masslimit* keywords cannot both be
used in the same *reaxff/species* fix.  The *delete_rate_limit*
keyword can enforce an upper limit on the overall rate of molecule
deletion.  The number of deletion occurrences is limited to Nlimit
within an interval of Nsteps timesteps.   Nlimit can be specified with
an equal-style :doc:`variable <variable>`.  When using the
*delete_rate_limit* keyword, no deletions are permitted to occur
within the first Nsteps timesteps of the first run (after reading a
either a data or restart file).

----------

The *Nevery*, *Nrepeat*, and *Nfreq* arguments specify on what
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

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.

This fix computes both a global vector of length 2 and a per-atom vector,
either of which can be accessed by various :doc:`output commands <Howto_output>`.
The values in the global vector are "intensive".

The 2 values in the global vector are as follows:

* 1 = total number of molecules
* 2 = total number of distinct species

The per-atom vector stores the molecule ID for each atom as identified
by the fix.  If an atom is not in a molecule, its ID will be 0.
For atoms in the same molecule, the molecule ID for all of them
will be the same and will be equal to the smallest atom ID of
any atom in the molecule.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.
This fix is not invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The "fix reaxff/species" requires that :doc:`pair_style reaxff <pair_reaxff>` is used.
This fix is part of the REAXFF package.  It is only enabled if LAMMPS was built with that
package.  See the :doc:`Build package <Build_package>` page for more info.

To write gzipped species files, you must compile LAMMPS with the -DLAMMPS_GZIP option.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix reaxff/bonds <fix_reaxff_bonds>`

Default
"""""""

The default values for bond-order cutoffs are 0.3 for all I-J pairs.
The default element symbols are C, H, O, N.
Position files are not written by default.
