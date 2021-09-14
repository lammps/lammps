.. index:: fix acks2/reaxff
.. index:: fix acks2/reaxff/kk


fix acks2/reaxff command
====================

Accelerator Variants: *acks2/reaxff/kk*

Syntax
""""""

.. parsed-literal::

   fix ID group-ID acks2/reaxff Nevery cutlo cuthi tolerance params args

* ID, group-ID are documented in :doc:`fix <fix>` command
* acks2/reaxff = style name of this fix command
* Nevery = perform ACKS2 every this many steps
* cutlo,cuthi = lo and hi cutoff for Taper radius
* tolerance = precision to which charges will be equilibrated
* params = reaxff or a filename

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all acks2/reaxff 1 0.0 10.0 1.0e-6 reaxff
   fix 1 all acks2/reaxff 1 0.0 10.0 1.0e-6 param.acks2

Description
"""""""""""

Perform the atom-condensed Kohn-Sham DFT to second order (ACKS2) charge
equilibration method as described in :ref:`(Verstraelen) <Verstraelen>`.
ACKS2 impedes unphysical long-range charge transfer sometimes seen with
QEq (e.g. for dissociation of molecules), at increased computational cost.
It is typically used in conjunction with the ReaxFF force field model as
implemented in the :doc:`pair_style reaxff <pair_reaxff>` command, but
it can be used with any potential in LAMMPS, so long as it defines and
uses charges on each atom.  For more technical details about the
charge equilibration performed by fix acks2/reaxff, see the
:ref:`(O'Hearn) <O'Hearn>` paper.

The ACKS2 method minimizes the electrostatic energy of the system by
adjusting the partial charge on individual atoms based on interactions
with their neighbors. It requires some parameters
for each atom type.
If the *params* setting above is the word "reaxff", then these are
extracted from the :doc:`pair_style reaxff <pair_reaxff>` command and
the ReaxFF force field file it reads in.  If a file name is specified
for *params*\ , then the parameters are taken from the specified file
and the file must contain one line for each atom type.  The latter
form must be used when performing QeQ with a non-ReaxFF potential.
The lines should be formatted as follows:

.. parsed-literal::

   bond_softness
   itype chi eta gamma bcut

where the first line is the global parameter *bond_softness*. The remaining
1 to Ntypes lines include *itype*, the atom type from 1 to Ntypes, *chi*, the
electronegativity in eV, *eta*, the self-Coulomb
potential in eV, *gamma*, the valence orbital
exponent, and *bcut*, the bond cutoff distance.  Note that these 4 quantities are also in the ReaxFF
potential file, except that eta is defined here as twice the eta value
in the ReaxFF file. Note that unlike the rest of LAMMPS, the units
of this fix are hard-coded to be A, eV, and electronic charge.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  No global scalar or vector or per-atom
quantities are stored by this fix for access by various :doc:`output commands <Howto_output>`.  No parameter of this fix can be used
with the *start/stop* keywords of the :doc:`run <run>` command.

This fix is invoked during :doc:`energy minimization <minimize>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This fix is part of the REAXFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix does not correctly handle interactions
involving multiple periodic images of the same atom. Hence, it should not
be used for periodic cell dimensions less than 10 angstroms.

Related commands
""""""""""""""""

:doc:`pair_style reaxff <pair_reaxff>`, :doc:`fix qeq/reaxff <fix_qeq_reaxff>`

**Default:** none

----------

.. _O'Hearn:

**(O'Hearn)** O'Hearn, Alperen, Aktulga, SIAM J. Sci. Comput., 42(1), C1–C22 (2020).

.. _Verstraelen:

**(Verstraelen)** Verstraelen, Ayers, Speybroeck, Waroquier, J. Chem. Phys. 138, 074108 (2013).
