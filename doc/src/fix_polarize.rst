.. index:: fix polarize/bem/gmres
.. index:: fix polarize/bem/icc
.. index:: fix polarize/functional

fix polarize/bem/gmres command
==============================

fix polarize/bem/icc command
============================

fix polarize/functional command
===============================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID style nevery tolerance ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *polarize/bem/gmres* or *polarize/bem/icc* or *polarize/functional*
* Nevery = this fixed is invoked every this many timesteps
* tolerance = the relative tolerance for the iterative solver to stop


Examples
""""""""

.. code-block:: LAMMPS

   fix 2 interface polarize/bem/gmres 5 0.0001
   fix 1 interface polarize/bem/icc 1 0.0001
   fix 3 interface polarize/functional 1 0.001


Used in input scripts:

   .. parsed-literal::

      examples/PACKAGES/dielectric/in.confined
      examples/PACKAGES/dielectric/in.nopbc

Description
"""""""""""

These fixes compute induced charges at the interface between two
impermeable media with different dielectric constants. The interfaces
need to be discretized into vertices, each representing a boundary element.
The vertices are treated as if they were regular atoms or particles.
:doc:`atom_style dielectric <atom_style>` should be used since it defines
the additional properties of each interface particle such as
interface normal vectors, element areas, and local dielectric mismatch.
These fixes also require the use of :doc:`pair_style <pair_style>` and
:doc:`kspace_style <kspace_style>` with the *dielectric* suffix.
At every time step, given a configuration of the physical charges in the system
(such as atoms and charged particles) these fixes compute and update
the charge of the interface particles. The interfaces are allowed to move
during the simulation with appropriate time integrators (for example,
with :doc:`fix_rigid <fix_rigid>`).

Consider an interface between two media: one with dielectric constant
of 78 (water), the other of 4 (silica). The interface is discretized
into 2000 boundary elements, each represented by an interface particle. Suppose that
each interface particle has a normal unit vector pointing from the silica medium to water.
The dielectric difference along the normal vector is then 78 - 4 = 74,
the mean dielectric value is (78 + 4) / 2 = 41. Each boundary element
also has its area and the local mean curvature (which is used by these fixes
for computing a correction term in the local electric field).
To model charged interfaces, the interface particle will have a non-zero charge value,
coming from its area and surface charge density.

For non-interface particles such as atoms and charged particles,
the interface normal vectors, element area, and dielectric mismatch are
irrelevant. Their local dielectric value is used to rescale their actual charge
when computing the Coulombic interactions. For instance, for a cation carrying
a charge of +2 (in charge unit) in an implicit solvent with dielectric constant of 40
would have actual charge of +2, and a local dielectric constant value of 40.
It is assumed that the particles cannot pass through the interface during the simulation
so that its local dielectric constant value does not change.

There are some example scripts for using these fixes
with LAMMPS in the ``examples/PACKAGES/dielectric`` directory. The README file
therein contains specific details on the system setup. Note that the example data files
show the additional fields (columns) needed for :doc:`atom_style dielectric <atom_style>`
beyond the conventional fields *id*, *mol*, *type*, *q*, *x*, *y*, and *z*.

----------

For fix *polarize/bem/gmres* and fix *polarize/bem/icc* the induced
charges of the atoms in the specified group, which are the vertices on
the interface, are computed using the equation:

.. math::

  \sigma_b(\mathbf{s}) = \dfrac{1 - \bar{\epsilon}}{\bar{\epsilon}}
     \sigma_f(\mathbf{s}) - \epsilon_0 \dfrac{\Delta \epsilon}{\bar{\epsilon}}
     \mathbf{E}(\mathbf{s}) \cdot \mathbf{n}(\mathbf{s})

* :math:`\sigma_b` is the induced charge density at the interface vertex :math:`\mathbf{s}`.
* :math:`\bar{\epsilon}` is the mean dielectric constant at the interface vertex: :math:`\bar{\epsilon} = (\epsilon_1 + \epsilon_2)/2`.
* :math:`\Delta \epsilon` is the dielectric constant difference at the interface vertex: :math:`\Delta \epsilon = \epsilon_1 - \epsilon_2`
* :math:`\sigma_f` is the free charge density at the interface vertex
* :math:`\mathbf{E}(\mathbf{s})` is the electrical field at the vertex
* :math:`\mathbf{n}(\mathbf{s})` is the unit normal vector at the vertex pointing from medium with :math:`\epsilon_2` to that with :math:`\epsilon_1`

Fix *polarize/bem/gmres* employs the Generalized Minimum Residual (GMRES)
as described in :ref:`(Barros) <Barros>` to solve :math:`\sigma_b`.

Fix *polarize/bem/icc* employs the successive over-relaxation algorithm
as described in :ref:`(Tyagi) <Tyagi>` to solve :math:`\sigma_b`.

The iterative solvers would terminate either when the maximum relative change
in the induced charges in consecutive iterations is below the set tolerance,
or when the number of iterations reaches *iter_max* (see below).

Fix *polarize/functional* employs the energy functional variation approach
as described in :ref:`(Jadhao) <Jadhao>` to solve :math:`\sigma_b`.


More details on the implementation of these fixes and their recommended use
are described in :ref:`(NguyenTD) <NguyenTD>`.


Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` command provides certain options to
control the induced charge solver and the initial values of the interface elements:

  .. parsed-literal::
      *itr_max* arg
         arg = maximum number of iterations for convergence
      *dielectrics* ediff emean epsilon area charge
         ediff = dielectric difference
         emean = dielectric mean
         epsilon = local dielectric value
         aree = element area
         charge = real interface charge

*polarize/bem/gmres* or *polarize/bem/icc* compute a global 2-element vector
which can be accessed by various :doc:`output commands <Howto_output>`.
The first element is the number of iterations when the solver terminates
(of which the upper bound is set by *iter_max*). The second element is the RMS error.


Restrictions
""""""""""""

These fixes are part of the DIELECTRIC package.  It is only enabled
if LAMMPS was built with that package, which requires that also the
KSPACE package is installed.  See the :doc:`Build package
<Build_package>` page for more info.

Note that the *polarize/bem/gmres* and *polarize/bem/icc* fixes only support
:doc:`units <units>` *lj*, *real*, *metal*, *si* and *nano* at the moment.


Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix polarize <fix_polarize>`, :doc:`read_data <read_data>`,
:doc:`pair_style lj/cut/coul/long/dielectric <pair_dielectric>`,
:doc:`kspace_style pppm/dielectric <kspace_style>`,
:doc:`compute efield/atom <compute_efield_atom>`

Default
"""""""

*iter_max* = 20

----------

.. _Barros:

**(Barros)** Barros, Sinkovits, Luijten, J. Chem. Phys, 140, 064903 (2014)

.. _Tyagi:

**(Tyagi)** Tyagi, Suzen, Sega, Barbosa, Kantorovich, Holm, J Chem Phys, 132, 154112 (2010)

.. _Jadhao:

**(Jadhao)** Jadhao, Solis, Olvera de la Cruz, J Chem Phys, 138, 054119 (2013)

.. _NguyenTD:

**(NguyenTD)** Nguyen, Li, Bagchi, Solis, Olvera de la Cruz, Comput Phys Commun 241, 80-19 (2019)

