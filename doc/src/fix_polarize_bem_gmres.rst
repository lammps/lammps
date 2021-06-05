.. index:: fix polarize/bem/gmres

fix polarize/bem/gmres command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID polarize/bem/gmres nevery tolerance ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* polarize/bem/gmres = style name of this fix command
* Nevery = this fixed is invoked every this many timesteps
* tolerance = the tolerance for the iterative solver to stop


Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all polarize/bem/gmres 5 0.0001

Description
"""""""""""

The fix polarize/bem/gmres computes the induced charges
at the interface between two impermeable media with
different dielectric constants.

There are some example scripts for using this fix
with LAMMPS in the examples/USER/dielectric directory.

----------

The induced charges of the atoms in the specified group, which are
the vertices on the interface, are computed using the equation:

..math::

  \sigma_b(\mathbf{s}) = \dfrac{1 - \bar{\epsilon}}{\bar{\epsilon}}
     \sigma_f(\mathbf{s}) - \epsilon_0 \dfrac{\Delta \epsilon}{\bar{\epsilon}}
     \mathbf{E}(\mathbf{s}) \cdot \mathbf{n}(\mathbf{s})

* :math:`\sigma_b` is the induced charge density at the interface vertex
:math:`\mathbf{s}`.
* :math:`\bar{\epsilon}` is the mean dielectric constant at the interface vertex:
:math:`\bar{\epsilon} = (\epsilon_1 + \epsilon_2)/2`.
* :math:`\Delta \epsilon` is the dielectric constant difference at the interface vertex:
:math:`\Delta \epsilon = \epsilon_1 - \epsilon_2`
* :math:`\sigma_f` is the free charge density at the interface vertex
* :math:`\mathbf{E}(\mathbf{s})` is the electrical field at the vertex
* :math:`\mathbf{n}(\mathbf{s})` is the unit normal vector at the vertex
pointing from medium with :math:`\epsilon_2` to that with :math:`\epsilon_1`


The fix polarize/bem/gmres employs the Generalized Minimum Residual (GMRES)
as described in :ref:`(Barros) <Barros>` to solve for :math:`\sigma_b`.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`fix_modify <fix_modify>`

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the Colvars
energy mentioned above.  The scalar value calculated by this fix is
"extensive".

Restrictions
""""""""""""

This fix is part of the USER-DIELECTRIC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

There can only be one colvars fix active at a time. Since the interface
communicates only the minimum amount of information and colvars module
itself can handle an arbitrary number of collective variables, this is
not a limitation of functionality.

Related commands
""""""""""""""""

:doc:`fix polarize/bem/icc <fix_polarize_bem_icc>`, :doc:`fix polarize/functional <fix_polarize_functional>`

Default
"""""""

None.

----------

.. _Barros:

**(Barros)** Barros, Sinkovits, Luijten, J. Chem. Phys, 140, 064903 (2014)


.. _NguyenTD:

**(NguyenTD)** Nguyen, Li, Bagchi, Solis, Olvera de la Cruz, Comput Phys Commun 241, 80-19 (2019)

