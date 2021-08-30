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
* tolerance = the tolerance for the iterative solver to stop


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
impermeable media with different dielectric constants.

There are some example scripts for using this fix
with LAMMPS in the examples/PACKAGES/dielectric directory.

----------

For fix *polarize/bem/gmres* and fix *polarize/bem/icc* the induced
charges of the atoms in the specified group, which are the vertices on
the interface, are computed using the equation:

..math::

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

Fix *polarize/functional* ...

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

...

Restrictions
""""""""""""

These fixes are part of the DIELECTRIC package.  It is only enabled
if LAMMPS was built with that package, which requires that also the
KSPACE package is installed.  See the :doc:`Build package
<Build_package>` page for more info.


Related commands
""""""""""""""""

:doc:`compute efield/atom <compute_efield_atom>`

Default
"""""""

None.

----------

.. _Barros:

**(Barros)** Barros, Sinkovits, Luijten, J. Chem. Phys, 140, 064903 (2014)

.. _Tyagi:

**(Tyagi)** Tyagi, Suzen, Sega, Barbosa, Kantorovich, Holm, J Chem Phys, 132, 154112 (2010)

.. _Jadhao:

**(Jadhao)** Jadhao, Solis, Olvera de la Cruz, J Chem Phys, 138, 054119 (2013)

.. _NguyenTD:

**(NguyenTD)** Nguyen, Li, Bagchi, Solis, Olvera de la Cruz, Comput Phys Commun 241, 80-19 (2019)

