Citing LAMMPS
-------------

Core Algorithms
^^^^^^^^^^^^^^^

Since LAMMPS is a community project, there is not a single one
publication or reference that describes **all** of LAMMPS.
The canonical publication that describes the foundation, that is
the basic spatial decomposition approach, the neighbor finding,
and basic communications algorithms used in LAMMPS is:

 `S. Plimpton, Fast Parallel Algorithms for Short-Range Molecular Dynamics, J Comp Phys, 117, 1-19 (1995). <http://www.sandia.gov/~sjplimp/papers/jcompphys95.pdf>`_

So any project using LAMMPS (or a derivative application using LAMMPS as
a simulation engine) should cite this paper. A new publication
describing the developments and improvements of LAMMPS in the 25 years
since then is currently in preparation.


DOI for the LAMMPS code
^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS developers use the `Zenodo service at CERN
<https://zenodo.org/>`_ to create digital object identifies (DOI) for
stable releases of the LAMMPS code. There are two types of DOIs for the
LAMMPS source code: the canonical DOI for **all** versions of LAMMPS,
which will always point to the **latest** stable release version is:

- DOI: `10.5281/zenodo.3726416 <https://dx.doi.org/10.5281/zenodo.3726416>`_

In addition there are DOIs for individual stable releases. Currently there are:

- 3 March 2020 version: `DOI:10.5281/zenodo.3726417 <https://dx.doi.org/10.5281/zenodo.3726417>`_
- 29 October 2020 version: `DOI:10.5281/zenodo.4157471 <https://dx.doi.org/10.5281/zenodo.4157471>`_


Home page
^^^^^^^^^

The LAMMPS website at `https://www.lammps.org/
<https://www.lammps.org>`_ is the canonical location for information
about LAMMPS and its features.

Citing contributions
^^^^^^^^^^^^^^^^^^^^

LAMMPS has many features and that use either previously published
methods and algorithms or novel features.  It also includes potential
parameter filed for specific models.  Where available, a reminder about
references for optional features used in a specific run is printed to
the screen and log file.  Style and output location can be selected with
the :ref:`-cite command-line switch <cite>`.  Additional references are
given in the documentation of the :doc:`corresponding commands
<Commands_all>` or in the :doc:`Howto tutorials <Howto>`.
