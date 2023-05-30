Citing LAMMPS
-------------

Core Algorithms
^^^^^^^^^^^^^^^

The paper mentioned below is the best overview of LAMMPS, but there are
also publications describing particular models or algorithms implemented
in LAMMPS or complementary software that is has interfaces to.  Please
see below for how to cite contributions to LAMMPS.

.. _lammps_paper:

The latest canonical publication that describes the basic features, the
source code design, the program structure, the spatial decomposition
approach, the neighbor finding, basic communications algorithms, and how
users and developers have contributed to LAMMPS is:

  `LAMMPS - A flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales, Comp. Phys. Comm. 271, 108171 (2022) <https://doi.org/10.1016/j.cpc.2021.108171>`_

So a project using LAMMPS or a derivative application that uses LAMMPS
as a simulation engine should cite this paper.  The paper is expected to
be published in its final form under the same DOI in the first half
of 2022.  Please also give the URL of the LAMMPS website in your paper,
namely https://www.lammps.org.

The original publication describing the parallel algorithms used in the
initial versions of LAMMPS is:

  `S. Plimpton, Fast Parallel Algorithms for Short-Range Molecular Dynamics, J Comp Phys, 117, 1-19 (1995). <https://doi.org/10.1006/jcph.1995.1039>`_


DOI for the LAMMPS source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The LAMMPS developers use the `Zenodo service at CERN <https://zenodo.org/>`_
to create digital object identifiers (DOI) for stable releases of the
LAMMPS source code.  There are two types of DOIs for the LAMMPS source code.

The canonical DOI for **all** versions of LAMMPS, which will always
point to the **latest** stable release version, is:

- DOI: `10.5281/zenodo.3726416 <https://dx.doi.org/10.5281/zenodo.3726416>`_

In addition there are DOIs generated for individual stable releases:

- 3 March 2020 version: `DOI:10.5281/zenodo.3726417 <https://dx.doi.org/10.5281/zenodo.3726417>`_
- 29 October 2020 version: `DOI:10.5281/zenodo.4157471 <https://dx.doi.org/10.5281/zenodo.4157471>`_
- 29 September 2021 version: `DOI:10.5281/zenodo.6386596 <https://dx.doi.org/10.5281/zenodo.6386596>`_

Home page
^^^^^^^^^

The LAMMPS website at `https://www.lammps.org/
<https://www.lammps.org>`_ is the canonical location for information
about LAMMPS and its features.

Citing contributions
^^^^^^^^^^^^^^^^^^^^

LAMMPS has many features that use either previously published methods
and algorithms or novel features.  It also includes potential parameter
files for specific models.  Where available, a reminder about references
for optional features used in a specific run is printed to the screen
and log file.  Style and output location can be selected with the
:ref:`-cite command-line switch <cite>`.  Additional references are
given in the documentation of the :doc:`corresponding commands
<Commands_all>` or in the :doc:`Howto tutorials <Howto>`.  Please make
certain, that you provide the proper acknowledgments and citations in
any published works using LAMMPS.
