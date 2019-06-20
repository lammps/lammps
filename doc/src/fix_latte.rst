.. index:: fix latte

fix latte command
=================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID latte peID

* ID, group-ID are documented in :doc:`fix <fix>` command
* latte = style name of this fix command
* peID = NULL or ID of compute used to calculate per-atom energy

Examples
""""""""


.. parsed-literal::

   fix dftb all latte NULL

Description
"""""""""""

This fix style is a wrapper on the self-consistent charge transfer
density functional based tight binding (DFTB) code LATTE. If you
download and build LATTE, it can be called as a library by LAMMPS via
this fix to run dynamics or perform energy minimization using DFTB
forces and energies computed by LATTE.

LATTE is principally developed and supported by Marc Cawkwell and
co-workers at Los Alamos National Laboratory (LANL).  See the full
list of contributors in the src/LATTE/README file.

To use this fix, the LATTE program needs to be compiled as a library
and linked with LAMMPS.  LATTE can be downloaded (or cloned) from
`https://github.com/lanl/LATTE <https://github.com/lanl/LATTE>`_.
Instructions on how to download and build LATTE on your system can be
found in the lib/latte/README.  Note that you can also use the "make
lib-latte" command from the LAMMPS src directory to automate this
process.

Once LAMMPS is built with the LATTE package, you can run the example
input scripts for molecular dynamics or energy minimization that are
found in examples/latte.

A step-by-step tutorial can be followed at: `LAMMPS-LATTE tutorial <https://github.com/lanl/LATTE/wiki/Using-LATTE-through-LAMMPS>`_

The *peID* argument is not yet supported by fix latte, so it must be
specified as NULL.  Eventually it will be used to enable LAMMPS to
calculate a Coulomb potential as an alternative to LATTE performing
the calculation.


----------


LATTE is a code for performing self-consistent charge transfer
tight-binding (SC-TB) calculations of total energies and the forces
acting on atoms in molecules and solids. This tight-binding method is
becoming more and more popular and widely used in chemistry,
biochemistry, material science, etc.

The SC-TB formalism is derived from an expansion of the Kohn-Sham
density functional to second order in charge fluctuations about a
reference charge of overlapping atom-centered densities and bond
integrals are parameterized using a Slater-Koster tight-binding
approach. This procedure, which usually is referred to as the DFTB
method has been described in detail by (:ref:`Elstner <Elstner>`) and
(:ref:`Finnis <Finnis2>`) and coworkers.

The work of the LATTE developers follows that of Elstner closely with
respect to the physical model.  However, the development of LATTE is
geared principally toward large-scale, long duration, microcanonical
quantum-based Born-Oppenheimer molecular dynamics (QMD) simulations.
One of the main bottlenecks of an electronic structure calculation is
the solution of the generalized eigenvalue problem which scales with
the cube of the system size O(N\^3).

The Theoretical and Computer sciences divisions at Los Alamos National
Laboratory have accumulated large experience addressing this issue by
calculating the density matrix directly instead of using
diagonalization. We typically use a recursive sparse Fermi-operator
expansion using second-order spectral projection functions
(SP2-algorithm), which was introduced by Niklasson in 2002
(:ref:`Niklasson2002 <Niklasson2002>`), (:ref:`Rubensson <Rubensson>`),
(:ref:`Mniszewski <Mniszewski>`).  When the matrices involved in the
recursive expansion are sufficiently sparse, the calculation of the
density matrix scales linearly as a function of the system size O(N).

Another important feature is the extended Lagrangian framework for
Born-Oppenheimer molecular dynamics (XL-BOMD)
(:ref:`Niklasson2008 <Niklasson2008>`) (:ref:`Niklasson2014 <Niklasson2014>`),
(:ref:`Niklasson2017 <Niklasson2017>`) that allows for a drastic reduction
or even a complete removal of the iterative self-consistent field
optimization.  Often only a single density matrix calculation per
molecular dynamics time step is required, yet total energy stability
is well maintained.  The SP2 and XL-BOMD techniques enables stable
linear scaling MD simulations with a very small computational
overhead.  This opens a number of opportunities in many different
areas of chemistry and materials science, as we now can simulate
larger system sizes and longer time scales
(:ref:`Cawkwell2012 <Cawkwell2012>`), (:ref:`Negre2016 <Negre2016>`).


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the potential energy computed by LATTE to the system's
potential energy as part of :doc:`thermodynamic output <thermo_style>`.

The :doc:`fix\_modify <fix_modify>` *virial* option is supported by this
fix to add the LATTE DFTB contribution to the system's virial as part
of :doc:`thermodynamic output <thermo_style>`.  The default is *virial
yes*

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy discussed above.  The scalar value calculated by this fix is
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The DFTB forces computed by LATTE via this fix are imposed during an
energy minimization, invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the potential energy associated with the DFTB
   forces to be included in the total potential energy of the system (the
   quantity being minimized), you MUST enable the
   :doc:`fix\_modify <fix_modify>` *energy* option for this fix.

Restrictions
""""""""""""


This fix is part of the LATTE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

You must use metal units, as set by the :doc:`units <units>` command to
use this fix.

LATTE does not currently compute per-atom energy or per-atom virial
contributions.  So they will not show up as part of the calculations
performed by the :doc:`compute pe/atom <compute_pe_atom>` or :doc:`compute stress/atom <compute_stress_atom>` commands.

Currently, LAMMPS must be run in serial or as a single MPI task, to
use this fix.  This is typically not a bottleneck, since LATTE will be
doing 99% or more of the work to compute quantum-accurate forces.

.. note::

   NEB calculations can be done using this fix using multiple
   replicas and running LAMMPS in parallel.  However, each replica must
   be run on a single MPI task.  For details, see the :doc:`neb <neb>`
   command doc page and the :doc:`-partition command-line switch <Run_options>`

**Related commands:** none

**Default:** none


----------


.. _Elstner:



**(Elstner)** M. Elstner, D. Poresag, G. Jungnickel, J. Elsner,
M. Haugk, T. Frauenheim, S. Suhai, and G. Seifert, Phys. Rev. B, 58,
7260 (1998).

.. _Elstner1:



**(Elstner)** M. Elstner, D. Poresag, G. Jungnickel, J. Elsner,
M. Haugk, T. Frauenheim, S. Suhai, and G. Seifert, Phys. Rev. B, 58,
7260 (1998).

.. _Finnis2:



**(Finnis)** M. W. Finnis, A. T. Paxton, M. Methfessel, and M. van
Schilfgarde, Phys. Rev. Lett., 81, 5149 (1998).

.. _Mniszewski:



**(Mniszewski)** S. M. Mniszewski, M. J. Cawkwell, M. E. Wall,
J. Mohd-Yusof, N. Bock, T. C.  Germann, and A. M. N. Niklasson,
J. Chem. Theory Comput., 11, 4644 (2015).

.. _Niklasson2002:



**(Niklasson2002)** A. M. N. Niklasson, Phys. Rev. B, 66, 155115 (2002).

.. _Rubensson:



**(Rubensson)** E. H. Rubensson, A. M. N. Niklasson, SIAM
J. Sci. Comput. 36 (2), 147-170, (2014).

.. _Niklasson2008:



**(Niklasson2008)** A. M. N. Niklasson, Phys. Rev. Lett., 100, 123004
(2008).

.. _Niklasson2014:



**(Niklasson2014)** A. M. N. Niklasson and M. Cawkwell, J. Chem. Phys.,
141, 164123, (2014).

.. _Niklasson2017:



**(Niklasson2017)** A. M. N. Niklasson, J. Chem. Phys., 147, 054103 (2017).

.. _Cawkwell2012:



**(Cawkwell2012)** A. M. N. Niklasson, M. J. Cawkwell, Phys. Rev. B, 86
(17), 174308 (2012).

.. _Negre2016:



**(Negre2016)** C. F. A. Negre, S. M. Mniszewski, M. J. Cawkwell,
N. Bock, M. E. Wall, and A. M. N. Niklasson, J. Chem. Theory Comp.,
12, 3063 (2016).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
