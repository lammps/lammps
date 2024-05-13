.. index:: pair_style lj/cut/tip4p/cut
.. index:: pair_style lj/cut/tip4p/cut/omp
.. index:: pair_style lj/cut/tip4p/long
.. index:: pair_style lj/cut/tip4p/long/gpu
.. index:: pair_style lj/cut/tip4p/long/omp
.. index:: pair_style lj/cut/tip4p/long/opt

pair_style lj/cut/tip4p/cut command
===================================

Accelerator Variants: *lj/cut/tip4p/cut/omp*

pair_style lj/cut/tip4p/long command
====================================

Accelerator Variants: *lj/cut/tip4p/long/gpu*, *lj/cut/tip4p/long/omp*, *lj/cut/tip4p/long/opt*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut/tip4p/cut* or *lj/cut/tip4p/long*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/cut/tip4p/cut* args = otype htype btype atype qdist cutoff (cutoff2)
       otype,htype = atom types for TIP4P O and H
       btype,atype = bond and angle types for TIP4P waters
       qdist = distance from O atom to massless charge (distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/tip4p/long* args = otype htype btype atype qdist cutoff (cutoff2)
       otype,htype = atom types for TIP4P O and H
       btype,atype = bond and angle types for TIP4P waters
       qdist = distance from O atom to massless charge (distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/tip4p/cut 1 2 7 8 0.15 12.0
   pair_style lj/cut/tip4p/cut 1 2 7 8 0.15 12.0 10.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/tip4p/long 1 2 7 8 0.15 12.0
   pair_style lj/cut/tip4p/long 1 2 7 8 0.15 12.0 10.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

Description
"""""""""""

The *lj/cut/tip4p* styles implement the TIP4P water model of
:ref:`(Jorgensen) <Jorgensen2>` and similar models, which introduce a
massless site M located a short distance away from the oxygen atom along
the bisector of the HOH angle.  The atomic types of the oxygen and
hydrogen atoms, the bond and angle types for OH and HOH interactions,
and the distance to the massless charge site are specified as pair_style
arguments and are used to identify the TIP4P-like molecules and
determine the position of the M site from the positions of the hydrogen
and oxygen atoms of the water molecules.  The M site location is used
for all Coulomb interactions instead of the oxygen atom location, also
with all other atom types, while the location of the oxygen atom is used
for the Lennard-Jones interactions.  Style *lj/cut/tip4p/cut* uses a
cutoff for Coulomb interactions; style *lj/cut/tip4p/long* is for use
with a long-range Coulombic solver (Ewald or PPPM).

.. note::

   For each TIP4P water molecule in your system, the atom IDs for
   the O and 2 H atoms must be consecutive, with the O atom first.  This
   is to enable LAMMPS to "find" the 2 H atoms associated with each O
   atom.  For example, if the atom ID of an O atom in a TIP4P water
   molecule is 500, then its 2 H atoms must have IDs 501 and 502.

See the :doc:`Howto tip4p <Howto_tip4p>` page for more information
on how to use the TIP4P pair styles and lists of parameters to set.
Note that the neighbor list cutoff for Coulomb interactions is
effectively extended by a distance 2\*qdist when using the TIP4P pair
style, to account for the offset distance of the fictitious charges on
O atoms in water molecules.  Thus it is typically best in an
efficiency sense to use a LJ cutoff >= Coulombic cutoff + 2\*qdist, to
shrink the size of the neighbor list.  This leads to slightly larger
cost for the long-range calculation, so you can test the trade-off for
your model.

The *lj/cut/tip4p* styles compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff.

They add Coulombic pairwise interactions given by

.. math::

   E = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c

where :math:`C` is an energy-conversion constant, :math:`q_i` and :math:`q_j`
are the charges on the 2 atoms, and :math:`\epsilon` is the dielectric
constant which can be set by the :doc:`dielectric <dielectric>` command.
If one cutoff is specified in the pair_style command, it is used for
both the LJ and Coulombic terms.  If two cutoffs are specified, they are
used as cutoffs for the LJ and Coulombic terms respectively.

Style *lj/cut/tip4p/long* compute the same
Coulombic interactions as style *lj/cut/tip4p/cut* except that an
additional damping factor is applied to the Coulombic term so it can
be used in conjunction with the :doc:`kspace_style <kspace_style>`
command and its *ewald* or *pppm* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.


Coefficients
""""""""""""

For all of the *lj/cut* pair styles, the following coefficients must
be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands, or by mixing as
described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* LJ cutoff (distance units)

Note that :math:`\sigma` is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at :math:`2^{\frac{1}{6}} \sigma`.

The last coefficient is optional.  If not specified, the global
LJ cutoff specified in the pair_style command is used.

For *lj/cut/tip4p/cut* and *lj/cut/tip4p/long* only the LJ cutoff can be
specified since a Coulombic cutoff cannot be specified for an individual I,J
type pair. All type pairs use the same global Coulombic cutoff specified in
the pair_style command.

.. warning::

   Because of how these pair styles implement the coulomb interactions
   by implicitly defining a fourth site for the negative charge
   of the TIP4P and similar water models, special care must be taken
   when using these pair styles with other computations that also use
   charges.  Unless they are specially set up to also handle the implicit
   definition of the 4th site, results are likely incorrect.  Example:
   :doc:`compute dipole/chunk <compute_dipole_chunk>`.  For the same
   reason, when using one of these pair styles with
   :doc:`pair_style  hybrid <pair_hybrid>`, **all** coulomb interactions
   should be handled by a single sub-style with TIP4P support. All other
   instances and styles will "see" the M point charges at the position
   of the Oxygen atom and thus compute incorrect forces and energies.
   LAMMPS will print a warning when it detects one of these issues.

----------

A version of these styles with a soft core, *lj/cut/tip4p/long/soft*, suitable
for use in free energy calculations, is part of the FEP package and
is documented with the :doc:`pair_style */soft <pair_fep_soft>`
styles.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair_modify" command
for details.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

The *lj/cut/coul/long* and *lj/cut/tip4p/long* pair styles support the
:doc:`pair_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` tail option for adding a long-range
tail correction to the energy and pressure for the Lennard-Jones
portion of the pair interaction.

All of the *lj/cut* pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

The *lj/cut* and *lj/cut/coul/long* pair styles support the use of the
*inner*, *middle*, and *outer* keywords of the :doc:`run_style respa <run_style>` command, meaning the pairwise forces can be
partitioned by distance at different levels of the rRESPA hierarchy.
The other styles only support the *pair* keyword of run_style respa.
See the :doc:`run_style <run_style>` command for details.

----------

Restrictions
""""""""""""

The *lj/cut/tip4p/long* styles are part of the
KSPACE package. The *lj/cut/tip4p/cut* style is part of the MOLECULE
package. These styles are only enabled if LAMMPS was built with those
packages.  See the :doc:`Build package <Build_package>` page for
more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Jorgensen2:

**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).
