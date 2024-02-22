.. index:: fix sgcmc

fix sgcmc command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID sgcmc every_nsteps swap_fraction temperature deltamu ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* sgcmc = style name of this fix command
* every_nsteps = number of MD steps between MC cycles
* swap_fraction = fraction of a full MC cycle carried out at each call (a value of 1.0 will perform as many trial moves as there are atoms)
* temperature = temperature that enters Boltzmann factor in Metropolis criterion (usually the same as MD temperature)
* deltamu = chemical potential difference(s) (`N-1` values must be provided, with `N` being the number of elements)
* Zero or more keyword/value pairs may be appended to fix definition line:

  .. parsed-literal::

     keyword = *variance* or *randseed* or *window_moves* or *window_size*
       *variance* kappa conc1 [conc2] ... [concN]
         kappa = variance constraint parameter
         conc1,conc2,...  = target concentration(s) in the range 0.0-1.0 (*N-1* values must be provided, with *N* being the number of elements)
       *randseed* N
         N = seed for pseudo random number generator
       *window_moves* N
         N = number of times sampling window is moved during one MC cycle
       *window_size* frac
         frac = size of sampling window (must be between 0.5 and 1.0)


Examples
""""""""

.. code-block:: LAMMPS

   fix mc all sgcmc 50 0.1 400.0 -0.55
   fix vc all sgcmc 20 0.2 700.0 -0.7 randseed 324234 variance 2000.0 0.05
   fix 2  all sgcmc 20 0.1 700.0 -0.7 window_moves 20

Description
"""""""""""

.. versionadded:: 22Dec2022

This command allows to carry out parallel hybrid molecular
dynamics/Monte Carlo (MD/MC) simulations using the algorithms described
in :ref:`(Sadigh1) <Sadigh1>`.  Simulations can be carried out in either
the semi-grand canonical (SGC) or variance constrained semi-grand
canonical (VC-SGC) ensemble :ref:`(Sadigh2) <Sadigh2>`. Only atom type
swaps are performed by the SGCMC fix. Relaxations are accounted for by
the molecular dynamics integration steps.

This fix can be used with standard multi-element EAM potentials
(:doc:`pair styles eam/alloy or eam/fs <pair_eam>`)

The SGCMC fix can handle Finnis/Sinclair type EAM potentials where
:math:`\rho(r)` is atom-type specific, such that different elements can
contribute differently to the total electron density at an atomic site
depending on the identity of the element at that atomic site.

------------

If this fix is applied, the regular MD simulation will be interrupted in
defined intervals to carry out a fraction of a Monte Carlo (MC)
cycle. The interval is set using the parameter *every_nsteps* which
determines how many MD integrator steps are taken between subsequent
calls to the MC routine.

It is possible to carry out pure lattice MC simulations by setting
*every_nsteps* to 1 and not defining an integration fix such as NVE,
NPT etc.  In that case, the particles will not move and only the MC
routine will be called to perform atom type swaps.

The parameter *swap_fraction* determines how many MC trial steps are carried
out every time the MC routine is entered. It is measured in units of full MC
cycles where one full cycle, *swap_fraction=1*, corresponds to as many MC
trial steps as there are atoms.

------------

The parameter *temperature* specifies the temperature that is used
to evaluate the Metropolis acceptance criterion. While it usually
should be set to the same value as the MD temperature there are cases
when it can be useful to use two different values for at least part of
the simulation, e.g., to speed up equilibration at low temperatures.

------------

The parameter *deltamu* is used to set the chemical potential difference
in the SGC MC algorithm (see Eq. 16 in :ref:`Sadigh1 <Sadigh1>`). By
convention it is the difference of the chemical potentials of elements
`B`, `C` ..., with respect to element A. When the simulation includes
`N` elements, `N-1` values must be specified.

------------

The variance-constrained SGC MC algorithm is activated if the keyword
*variance* is used. In that case the fix parameter *deltamu* determines
the effective average constraint in the parallel VC-SGC MC algorithm
(parameter :math:`\delta\mu_0` in Eq. (20) of :ref:`Sadigh1
<Sadigh1>`). The parameter *kappa* specifies the variance constraint
(see Eqs. (20-21) in :ref:`Sadigh1 <Sadigh1>`).

The parameter *conc* sets the target concentration (parameter
:math:`c_0` in Eqs.  (20-21) of :ref:`Sadigh1 <Sadigh1>`). The atomic
concentrations refer to components `B`, `C` ..., with `A` being set
automatically. When the simulation includes `N` elements, `N-1`
concentration values must be specified.

------------

There are several technical parameters that can be set via optional flags.

*randseed* is expected to be a positive integer number and is used
to initialize the random number generator on each processor.

*window_size* controls the size of the sampling window in a parallel MC
simulation. The size has to lie between 0.5 and 1.0. Normally, this
parameter should be left unspecified which instructs the code to choose
the optimal window size automatically (see Sect. III.B and Figure 6 in
:ref:`Sadigh1 <Sadigh1>` for details).

The number of times the window is moved during a MC cycle is set using
the parameter *window_moves* (see Sect. III.B in :ref:`Sadigh1
<Sadigh1>` for details).

------------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to restart files.

The MC routine keeps track of the global concentration(s) as well as the
number of accepted and rejected trial swaps during each MC step. These
values are provided by the sgcmc fix in the form of a global vector that
can be accessed by various :doc:`output commands <Howto_output>`
components of the vector represent the following quantities:

* 1 = The absolute number of accepted trial swaps during the last MC step
* 2 = The absolute number of rejected trial swaps during the last MC step
* 3 = The current global concentration of species *A* (= number of atoms of type 1 / total number of atoms)
* 4 = The current global concentration of species *B* (= number of atoms of type 2 / total number of atoms)
* ...
* N+2: The current global concentration of species *X* (= number of atoms of type *N* / total number of atoms)

Restrictions
""""""""""""

This fix is part of the MC package. It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

This fix style requires an :doc:`atom style <atom_style>` with per atom
type masses.

At present the fix provides optimized subroutines for EAM type
potentials (see above) that calculate potential energy changes due to
*local* atom type swaps very efficiently.  Other potentials are
supported by using the generic potential functions. This, however, will
lead to exceedingly slow simulations since it implies that the
energy of the *entire* system is recomputed at each MC trial step.  If
other potentials are to be used it is strongly recommended to modify and
optimize the existing generic potential functions for this purpose.
Also, the generic energy calculation can not be used for parallel
execution i.e. it only works with a single MPI process.

------------

Default
"""""""

The optional parameters default to the following values:

* *randseed* = 324234
* *window_moves* = 8
* *window_size* = automatic

------------

.. _Sadigh1:

**(Sadigh1)** B. Sadigh, P. Erhart, A. Stukowski, A. Caro, E. Martinez, and L. Zepeda-Ruiz, Phys. Rev. B **85**, 184203 (2012)

.. _Sadigh2:

**(Sadigh2)** B. Sadigh and P. Erhart, Phys. Rev. B **86**, 134204 (2012)
