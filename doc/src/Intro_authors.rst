Authors of LAMMPS
-----------------

The current core LAMMPS developers are listed here (grouped by seniority
and sorted alphabetically by last name). You can email an individual
developer with code related questions for their area of expertise, or
send an email to all of them at this address: "developers at
lammps.org".  General questions about LAMMPS should be posted in the
`LAMMPS forum on MatSci <https://matsci.org/lammps/>`_.

.. raw:: latex

   \small

.. list-table::
   :widths: 17 15 25 43
   :header-rows: 1

   * - Name
     - Affiliation
     - Email
     - Areas of expertise
   * - `Axel Kohlmeyer <ak_>`_
     - Temple U
     - akohlmey at gmail.com
     - OpenMP, library interfaces, LAMMPS-GUI, GitHub, MatSci forum, code maintenance, testing, releases
   * - `Steve Plimpton <sjp_>`_
     - SNL (retired)
     - sjplimp at gmail.com
     - MD kernels, parallel algorithms & scalability, code structure and design
   * - `Aidan Thompson <at_>`_
     - SNL
     - athomps at sandia.gov
     - manybody potentials, machine learned potentials, materials science, statistical mechanics
   * -
     -
     -
     -
   * - `Richard Berger <rb_>`_
     - LANL
     - richard.berger at outlook.com
     - Python, HPC, DevOps
   * - `Germain Clavier <gc_>`_
     - U Caen
     - germain.clavier at unicaen.fr
     - organic molecules, polymers, mechanical properties, surfaces, integrators, coarse-graining
   * - Joel Clemmer
     - SNL
     - jtclemm at sandia.gov
     -  granular systems fluid/solid mechanics
   * - `Jacob R. Gissinger <jg_>`_
     - Stevens Institute of Technology
     - jgissing at stevens.edu
     - reactive molecular dynamics, macro-molecular systems, type labels
   * - James Goff
     - SNL
     - jmgoff at sandia.gov
     - machine learned potentials, QEq solvers, Python
   * - Megan McCarthy
     - SNL
     - megmcca at sandia.gov
     - alloys, micro-structure, machine learned potentials
   * - Stan Moore
     - SNL
     - stamoor at sandia.gov
     - Kokkos, KSpace solvers, ReaxFF
   * - `Trung Nguyen <tn_>`_
     - U Chicago
     - ndactrung at gmail.com
     - soft matter, GPU package

.. _rb:  https://rbberger.github.io/
.. _gc:  https://enthalpiste.fr/
.. _jg:  https://www.nanocipher.org/
.. _ak:  https://sites.google.com/site/akohlmey/
.. _tn:  https://sites.google.com/site/ndtrung8/
.. _at:  https://www2.sandia.gov/~athomps/
.. _sjp: https://sjplimp.github.io
.. _lws: https://www.lammps.org

.. raw:: latex

   \normalsize

Past developers include Paul Crozier and Mark Stevens, both at SNL,
and Ray Shan, now at Materials Design.

----------

The `Authors page <https://www.lammps.org/authors.html>`_ of the
`LAMMPS website <lws_>`_ has a comprehensive list of all the individuals
who have contributed code for a new feature or command or tool to
LAMMPS.

----------

The following folks deserve special recognition.  Many of the packages
they have written are unique for an MD code and LAMMPS would not be as
general-purpose as it is without their expertise and efforts.

* Metin Aktulga (MSU), REAXFF package for C/C++ version of ReaxFF
* Mike Brown (Intel), GPU and INTEL packages
* Colin Denniston (U Western Ontario), LATBOLTZ package
* Georg Ganzenmuller (EMI), MACHDYN and SPH packages
* Andres Jaramillo-Botero (Caltech), EFF package for electron force field
* Reese Jones (Sandia) and colleagues, ATC package for atom/continuum coupling
* Christoph Kloss (DCS Computing), LIGGGHTS code for granular materials, built on top of LAMMPS
* Rudra Mukherjee (JPL), POEMS package for articulated rigid body motion
* Trung Ngyuen (U Chicago), GPU, RIGID, BODY, and DIELECTRIC packages
* Mike Parks (Sandia), PERI package for Peridynamics
* Roy Pollock (LLNL), Ewald and PPPM solvers
* Julien Tranchida (CEA Cadarache), SPIN package
* Christian Trott (Sandia), CUDA and KOKKOS packages
* Ilya Valuev (JIHT), AWPMD package for wave packet MD
* Greg Wagner (Northwestern U), MEAM package for MEAM potential

----------

As discussed on the `History page <https://www.lammps.org/history.html>`_ of the website, LAMMPS
originated as a cooperative project between DOE labs and industrial
partners.  Folks involved in the design and testing of the original
version of LAMMPS were the following:

* John Carpenter (Mayo Clinic, formerly at Cray Research)
* Terry Stouch (Lexicon Pharmaceuticals, formerly at Bristol Myers Squibb)
* Steve Lustig (Dupont)
* Jim Belak and Roy Pollock (LLNL)
