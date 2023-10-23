.. index:: fix phonon

fix phonon command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID phonon N Noutput Nwait map_file prefix keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* phonon = style name of this fix command
* N = measure the Green's function every this many timesteps
* Noutput = output the dynamical matrix every this many measurements
* Nwait = wait this many timesteps before measuring
* map_file = *file* or *GAMMA*

  .. parsed-literal::

       *file* is the file that contains the mapping info between atom ID and the lattice indices.

  .. parsed-literal::

       *GAMMA* flags to treate the whole simulation box as a unit cell, so that the mapping
       info can be generated internally. In this case, dynamical matrix at only the gamma-point
       will/can be evaluated.

* prefix = prefix for output files
* one or none keyword/value pairs may be appended
* keyword = *sysdim* or *nasr*

  .. parsed-literal::

       *sysdim* value = d
         d = dimension of the system, usually the same as the MD model dimension
       *nasr* value = n
         n = number of iterations to enforce the acoustic sum rule

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all phonon 20 5000 200000 map.in LJ1D sysdim 1
   fix 1 all phonon 20 5000 200000 map.in EAM3D
   fix 1 all phonon 10 5000 500000 GAMMA EAM0D nasr 100

Description
"""""""""""

Calculate the dynamical matrix from molecular dynamics simulations
based on fluctuation-dissipation theory for a group of atoms.

Consider a crystal with :math:`N` unit cells in three dimensions labeled
:math:`l = (l_1, l_2, l_3)` where :math:`l_i` are integers.  Each unit cell is
defined by three linearly independent vectors :math:`\mathbf{a}_1`,
:math:`\mathbf{a}_2`, :math:`\mathbf{a}_3` forming a parallelepiped,
containing :math:`K` basis atoms labeled :math:`k`.

Based on fluctuation-dissipation theory, the force constant
coefficients of the system in reciprocal space are given by
(:ref:`Campana <Campana>` , :ref:`Kong <Kong>`)

.. math::

   \mathbf{\Phi}_{k\alpha,k^\prime \beta}(\mathbf{q}) = k_B T \mathbf{G}^{-1}_{k\alpha,k^\prime \beta}(\mathbf{q})

where :math:`\mathbf{G}` is the Green's functions coefficients given by

.. math::

   \mathbf{G}_{k\alpha,k^\prime \beta}(\mathbf{q}) = \left< \mathbf{u}_{k\alpha}(\mathbf{q}) \bullet \mathbf{u}_{k^\prime \beta}^*(\mathbf{q}) \right>

where :math:`\left< \ldots \right>` denotes the ensemble average, and

.. math::

   \mathbf{u}_{k\alpha}(\mathbf{q}) = \sum_l \mathbf{u}_{l k \alpha} \exp{(i\mathbf{qr}_l)}

is the :math:`\alpha` component of the atomic displacement for the :math:`k`
th atom in the unit cell in reciprocal space at :math:`\mathbf{q}`. In
practice, the Green's functions coefficients can also be measured
according to the following formula,

.. math::

   \mathbf{G}_{k\alpha,k^\prime \beta}(\mathbf{q}) =
   \left< \mathbf{R}_{k \alpha}(\mathbf{q}) \bullet \mathbf{R}^*_{k^\prime \beta}(\mathbf{q}) \right>
   - \left<\mathbf{R}\right>_{k \alpha}(\mathbf{q}) \bullet \left<\mathbf{R}\right>^*_{k^\prime \beta}(\mathbf{q})

where :math:`\mathbf{R}` is the instantaneous positions of atoms, and
:math:`\left<\mathbf{R}\right>` is the averaged atomic positions. It
gives essentially the same results as the displacement method and is
easier to implement in an MD code.

Once the force constant matrix is known, the dynamical matrix
:math:`\mathbf{D}` can then be obtained by

.. math::

   \mathbf{D}_{k\alpha, k^\prime\beta}(\mathbf{q}) =
   (m_k m_{k^\prime})^{-\frac{1}{2}} \mathbf{\Phi}_{k \alpha, k^\prime \beta}(\mathbf{q})

whose eigenvalues are exactly the phonon frequencies at :math:`\mathbf{q}`.

This fix uses positions of atoms in the specified group and calculates
two-point correlations.  To achieve this. the positions of the atoms
are examined every *Nevery* steps and are Fourier-transformed into
reciprocal space, where the averaging process and correlation
computation is then done.  After every *Noutput* measurements, the
matrix :math:`\mathbf{G}(\mathbf{q})` is calculated and inverted to
obtain the elastic stiffness coefficients.  The dynamical matrices are
then constructed and written to *prefix*\ .bin.timestep files in binary
format and to the file *prefix*\ .log for each wave-vector
:math:`\mathbf{q}`.

A detailed description of this method can be found in
(:ref:`Kong2011 <Kong2011>`).

The *sysdim* keyword is optional.  If specified with a value smaller
than the dimensionality of the LAMMPS simulation, its value is used
for the dynamical matrix calculation.  For example, using LAMMPS to
model a 2D or 3D system, the phonon dispersion of a 1D atomic chain
can be computed using *sysdim* = 1.

The *nasr* keyword is optional.  An iterative procedure is employed to
enforce the acoustic sum rule on :math:`\Phi` at :math:`\Gamma`, and the number
provided by keyword *nasr* gives the total number of iterations. For a
system whose unit cell has only one atom, *nasr* = 1 is sufficient;
for other systems, *nasr* = 10 is typically sufficient.

The *map_file* contains the mapping information between the lattice
indices and the atom IDs, which tells the code which atom sits at
which lattice point; the lattice indices start from 0. An auxiliary
code, `latgen <https://code.google.com/p/latgen>`_, can be employed to
generate the compatible map file for various crystals.

In case one simulates a non-periodic system, where the whole simulation
box is treated as a unit cell, one can set *map_file* as *GAMMA*, so
that the mapping info will be generated internally and a file is not
needed. In this case, the dynamical matrix at only the gamma-point
will/can be evaluated. Please keep in mind that fix-phonon is designed
for cyrstals, it will be inefficient and even degrade the performance
of LAMMPS in case the unit cell is too large.

The calculated dynamical matrix elements are written out in
:doc:`energy/distance\^2/mass <units>` units.  The coordinates for *q*
points in the log file is in the units of the basis vectors of the
corresponding reciprocal lattice.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix. You can use it to change the temperature compute from thermo_temp
to the one that reflects the true temperature of atoms in the group.

No global scalar or vector or per-atom quantities are stored by this
fix for access by various :doc:`output commands <Howto_output>`.

Instead, this fix outputs its initialization information (including
mapping information) and the calculated dynamical matrices to the file
*prefix*\ .log, with the specified *prefix*\ .  The dynamical matrices are
also written to files *prefix*\ .bin.timestep in binary format.  These
can be read by the post-processing tool in tools/phonon to compute the
phonon density of states and/or phonon dispersion curves.

No parameter of this fix can be used with the *start/stop* keywords
of the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix assumes a crystalline system with periodical lattice. The
temperature of the system should not exceed the melting temperature to
keep the system in its solid state.

This fix is part of the PHONON package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This fix requires LAMMPS be built with an FFT library.  See the :doc:`Build settings <Build_settings>` page for details.

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`,
:doc:`dynamical_matrix <dynamical_matrix>`

Default
"""""""

The option defaults are sysdim = the same dimension as specified by
the :doc:`dimension <dimension>` command, and nasr = 20.

----------

.. _Campana:

**(Campana)** C. Campana and
M. H. Muser, *Practical Green's function approach to the
simulation of elastic semi-infinite solids*, `Phys. Rev. B [74], 075420 (2006) <https://doi.org/10.1103/PhysRevB.74.075420>`_

.. _Kong:

**(Kong)** L.T. Kong, G. Bartels, C. Campana,
C. Denniston, and Martin H. Muser, *Implementation of Green's
function molecular dynamics: An extension to LAMMPS*, `Computer Physics Communications [180](6):1004-1010 (2009). <https://doi.org/10.1016/j.cpc.2008.12.035>`_

L.T. Kong, C. Denniston, and Martin H. Muser,
*An improved version of the Green's function molecular dynamics
method*, `Computer Physics Communications [182](2):540-541 (2011). <https://doi.org/10.1016/j.cpc.2010.10.006>`_

.. _Kong2011:

**(Kong2011)** L.T. Kong, *Phonon dispersion measured directly from
molecular dynamics simulations*, `Computer Physics Communications [182](10):2201-2207, (2011). <https://doi.org/10.1016/j.cpc.2011.04.019>`_
