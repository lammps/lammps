.. index:: compute gaussian/grid/local
.. index:: compute gaussian/grid/local/kk

compute gaussian/grid/local command
===================================

Accelerator Variants: *gaussian/grid/local/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID gaussian/grid/local grid nx ny nz rcutfac  R_1 R_2 ... sigma_1 sigma_2

* ID, group-ID are documented in :doc:`compute <compute>` command
* gaussian/grid/local = style name of this compute command
* *grid* values = nx, ny, nz, number of grid points in x, y, and z directions (positive integer)
* *rcutfac* = scale factor applied to all cutoff radii (positive real)
* *R_1, R_2,...* = list of cutoff radii, one for each type (distance units)
* *sigma_1, sigma_2,...* = Gaussian widths, one for each type (distance units)

Examples
""""""""

.. code-block:: LAMMPS

    compute mygrid all gaussian/grid/local grid 40 40 40 4.0 0.5 0.5 0.4 0.4

Description
"""""""""""

.. versionadded:: 4Feb2025

Define a computation that calculates a Gaussian representation of the ionic
structure. This representation is used for the efficient evaluation
of quantities related to the structure factor in a grid-based workflow,
such as the ML-DFT workflow MALA :ref:`(Ellis) <Ellis2021b>`, for which it was originally
implemented. Usage of the workflow is described in a separate publication :ref:`(Fiedler) <Fiedler2023>`.

For each LAMMPS type, a separate sum of Gaussians is calculated, using
a separate Gaussian broadening per type. The computation
is always performed on the numerical grid, no atom-based version of this
compute exists. The Gaussian representation can only be executed in a local
fashion, thus the output array only contains rows for grid points
that are local to the processor subdomain. The layout of the grid is the same
as for the see :doc:`sna/grid/local <compute_sna_atom>` command.

Namely, the array contains one row for each of the
local grid points, looping over the global index *ix* fastest,
then *iy*, and *iz* slowest.  Each row of the array contains
the global indexes *ix*, *iy*, and *iz* first, followed by the *x*, *y*,
and *z* coordinates of the grid point, followed by the values of the Gaussians
(one floating point number per type per grid point).

----------


.. include:: accel_styles.rst



----------

Output info
"""""""""""

Compute *gaussian/grid/local* evaluates a local array.
The array contains one row for each of the
local grid points, looping over the global index *ix* fastest,
then *iy*, and *iz* slowest.  The array contains math :math:`ntypes+6` columns,
where *ntypes* is the number of LAMMPS types. The first three columns are
the global indexes *ix*, *iy*, and *iz*, followed by the *x*, *y*,
and *z* coordinates of the grid point, followed by the *ntypes* columns
containing the values of the Gaussians for each type.

Restrictions
""""""""""""

These computes are part of the ML-SNAP package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute sna/grid/local <compute_sna_atom>`

----------

.. _Ellis2021b:

**(Ellis)** Ellis, Fiedler, Popoola, Modine, Stephens, Thompson, Cangi, Rajamanickam, `Phys. Rev. B, 104, 035120, (2021) <https://doi.org/10.1103/PhysRevB.104.035120>`_

.. _Fiedler2023:

**(Fiedler)** Fiedler, Modine, Schmerler, Vogel, Popoola, Thompson, Rajamanickam, and Cangi,
`npj Comp. Mater., 9, 115 (2023) <https://doi.org/10.1038/s41524-023-01070-z>`_

