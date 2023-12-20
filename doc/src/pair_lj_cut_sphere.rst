.. index:: pair_style lj/cut/sphere
.. index:: pair_style lj/cut/sphere/omp

pair_style lj/cut/sphere command
================================

Accelerator Variant: *lj/cut/sphere/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut/sphere*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/cut/sphere* args = cutoff ratio
       cutoff = global cutoff ratio for Lennard Jones interactions (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/sphere 2.5
   pair_coeff * * 1.0
   pair_coeff 1 1 1.1 2.8

Description
"""""""""""

.. versionadded:: 15Jun2023

The *lj/cut/sphere* style compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} -
       \left(\frac{\sigma_{ij}}{r}\right)^6 \right]  \qquad  r < r_c * \sigma_{ij}

:math:`r_c` is the cutoff ratio.

This is the same potential function used by the :doc:`lj/cut
<pair_lj>` pair style, but the :math:`\sigma_{ij}` parameter is not
set as a per-type parameter via the :doc:`pair_coeff command
<pair_coeff>`.  Instead it is calculated individually for each pair
using the per-atom diameter attribute of :doc:`atom_style sphere
<atom_style>` for the two atoms as :math:`\sigma_{i}` and
:math:`\sigma_{j}`; :math:`\sigma_{ij}` is then computed by the mixing
rule for pair coefficients as set by the :doc:`pair_modify mix
<pair_modify>` command (defaults to geometric mixing).  The cutoff is
not specified as a distance, but as ratio that is internally
multiplied by :math:`\sigma_{ij}` to obtain the actual cutoff for each
pair of atoms.

Note that :math:`\sigma_{ij}` is defined in the LJ formula above as the
zero-crossing distance for the potential, *not* as the energy minimum which
is at :math:`2^{\frac{1}{6}} \sigma_{ij}`.

.. admonition:: Notes on cutoffs, neighbor lists, and efficiency
   :class: note

   If your system is mildly polydisperse, meaning the ratio of the
   diameter of the largest particle to the smallest is less than 2.0,
   then the neighbor lists built by the code should be reasonably
   efficient.  Which means they will not contain too many particle
   pairs that do not interact.  However, if your system is highly
   polydisperse (ratio > 2.0), the neighbor list build and force
   computations may be inefficient.  There are two ways to try and
   speed up the simulations.

   The first is to assign atoms to different atom types so that atoms of
   each type are similar in size.  E.g. if particle diameters range from
   1 to 5 use 4 atom types, ensuring atoms of type 1 have diameters from
   1.0-2.0, type 2 from 2.0-3.0, etc.  This will reduce the number of
   non-interacting pairs in the neighbor lists and thus reduce the time
   spent on computing pairwise interactions.

   The second is to use the :doc:`neighbor multi <neighbor>` command
   which enabled a different algorithm for building neighbor lists.  This
   will also require that you assign multiple atom types according to
   diameters, but will in addition use a more efficient size-dependent
   strategy to construct the neighbor lists and thus reduce the time
   spent on building neighbor lists.

   Here are example input script commands using both ideas for a
   highly polydisperse system:

   .. code-block:: c++

      units           lj
      atom_style      sphere
      lattice         fcc 0.8442
      region          box block 0 10 0 10 0 10
      create_box      2 box
      create_atoms    1 box

      # create atoms with random diameters from bimodal distribution
      variable switch atom random(0.0,1.0,345634)
      variable diam atom (v_switch<0.75)*normal(0.4,0.075,325)+(v_switch>=0.7)*normal(1.2,0.2,453)
      set group all diameter v_diam

      # assign type 2 to atoms with diameter > 0.6
      variable large atom (2.0*radius)>0.6
      group large variable large
      set group large type 2

      pair_style      lj/cut/sphere 2.5
      pair_coeff      * * 1.0

      neighbor 0.3 multi

   Using multiple atom types speeds up the calculation for this example
   by more than a factor of 2, and using the multi-style neighbor list
   build causes an additional speedup of about 20 percent.

Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* LJ cutoff ratio (unitless) (optional)

The last coefficient is optional.  If not specified, the global LJ
cutoff ratio specified in the :doc:`pair_style command <pair_style>` is
used.

If a repulsive only LJ interaction is desired, the coefficient for the cutoff
ratio should be set to the minimum of the LJ potential using ``$(2.0^(1.0/6.0))``

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon coefficients and cutoff
ratio for the *lj/cut/sphere* pair style can be mixed.  The default mixing
style is *geometric*.  See the :doc:`pair_modify command <pair_modify>`
for details.

The *lj/cut/sphere* pair style supports the :doc:`pair_modify shift <pair_modify>`
option for the energy of the Lennard-Jones portion of the pair interaction.

The *lj/cut/sphere* pair style does *not* support the :doc:`pair_modify
<pair_modify>` tail option for adding a long-range tail corrections to
the energy and pressure.

The *lj/cut/sphere* pair style writes its information to :doc:`binary
restart files <restart>`, so pair_style and pair_coeff commands do not
need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does *not* support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *lj/cut/sphere* pair style is only enabled if LAMMPS was built with the
EXTRA-PAIR package.  See the :doc:`Build package <Build_package>` page
for more info.

The *lj/cut/sphere* pair style does not support the *sixthpower* mixing rule.

----------

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut <pair_lj>`
* :doc:`pair_style lj/expnd/sphere <pair_lj_expand_sphere>`

Default
"""""""

none
