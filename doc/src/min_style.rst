.. index:: min\_style

min\_style command
==================

Syntax
""""""


.. parsed-literal::

   min_style style

* style = *cg* or *hftn* or *sd* or *quickmin* or *fire* or *fire/old* or *spin* or *spin/cg* or *spin/lbfgs*

Examples
""""""""


.. parsed-literal::

   min_style cg
   min_style spin
   min_style fire

Description
"""""""""""

Choose a minimization algorithm to use when a :doc:`minimize
<minimize>` command is performed.

Style *cg* is the Polak-Ribiere version of the conjugate gradient (CG)
algorithm.  At each iteration the force gradient is combined with the
previous iteration information to compute a new search direction
perpendicular (conjugate) to the previous search direction.  The PR
variant affects how the direction is chosen and how the CG method is
restarted when it ceases to make progress.  The PR variant is thought
to be the most effective CG choice for most problems.

Style *hftn* is a Hessian-free truncated Newton algorithm.  At each
iteration a quadratic model of the energy potential is solved by a
conjugate gradient inner iteration.  The Hessian (second derivatives)
of the energy is not formed directly, but approximated in each
conjugate search direction by a finite difference directional
derivative.  When close to an energy minimum, the algorithm behaves
like a Newton method and exhibits a quadratic convergence rate to high
accuracy.  In most cases the behavior of *hftn* is similar to *cg*\ ,
but it offers an alternative if *cg* seems to perform poorly.  This
style is not affected by the :doc:`min_modify <min_modify>` command.

Style *sd* is a steepest descent algorithm.  At each iteration, the
search direction is set to the downhill direction corresponding to the
force vector (negative gradient of energy).  Typically, steepest
descent will not converge as quickly as CG, but may be more robust in
some situations.

Style *quickmin* is a damped dynamics method described in
:ref:`(Sheppard) <Sheppard>`, where the damping parameter is related
to the projection of the velocity vector along the current force
vector for each atom.  The velocity of each atom is initialized to 0.0
by this style, at the beginning of a minimization.

Style *fire* is a damped dynamics method described in :ref:`(Bitzek)
<Bitzek>`, which is similar to *quickmin* but adds a variable timestep
and alters the projection operation to maintain components of the
velocity non-parallel to the current force vector.  The velocity of
each atom is initialized to 0.0 by this style, at the beginning of a
minimization. This style correspond to an optimized version described
in :ref:`(Guenole) <Guenole>` that include different time integration
schemes and defaults parameters. The default parameters can be
modified with the command :doc:`min_modify <min_modify>`.


Style *fire/old* is the original implementation of *fire* in Lammps,
conserved for backward compatibility. The main differences regarding
the current version *fire* are: time integration by Explicit Euler
only, different sequence in maintaining velocity components non-parallel
to the current force vector and hard-coded minimization parameters.
A complete description of the differences between *fire/old* and *fire*
can be found in :ref:`(Guenole) <Guenole>` (where the current *fire*
in LAMMPS is called *fire2.0*). By using an appropriate set of
parameters, *fire* can behave similar to *fire/old*, as described
in the :doc:`min_modify <min_modify>` command.

Style *spin* is a damped spin dynamics with an adaptive timestep.

Style *spin/cg* uses an orthogonal spin optimization (OSO) combined to
a conjugate gradient (CG) approach to minimize spin configurations.

Style *spin/lbfgs* uses an orthogonal spin optimization (OSO) combined
to a limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) approach
to minimize spin configurations.

See the :doc:`min/spin <min_spin>` doc page for more information about
the *spin*\ , *spin/cg* and *spin/lbfgs* styles.

Either the *quickmin*\ , *fire* and *fire/old* styles are useful in the
context of nudged elastic band (NEB) calculations via the :doc:`neb
<neb>` command.

Either the *spin*\ , *spin/cg* and *spin/lbfgs* styles are useful in
the context of magnetic geodesic nudged elastic band (GNEB)
calculations via the :doc:`neb/spin <neb_spin>` command.

.. note::

   The damped dynamic minimizers use whatever timestep you have
   defined via the :doc:`timestep <timestep>` command.  Often they
   will converge more quickly if you use a timestep about 10x larger
   than you would normally use for dynamics simulations.
   For *fire*, the default timestep is recommended to be equal to
   the one you would normally use for dynamics simulations.

.. note::

   The *quickmin*\ , *fire*\ , *fire/old*\ , *hftn*\ , and *cg/kk* styles do not yet
   support the use of the :doc:`fix box/relax <fix_box_relax>` command
   or minimizations involving the electron radius in :doc:`eFF
   <pair_eff>` models.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>`
doc page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package
<Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix
command-line switch <Run_options>` when you invoke LAMMPS, or you can
use the :doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`min_modify <min_modify>`, :doc:`minimize <minimize>`, :doc:`neb <neb>`

Default
"""""""


.. parsed-literal::

   min_style cg


----------


.. _Sheppard:

**(Sheppard)** Sheppard, Terrell, Henkelman, J Chem Phys, 128, 134106
(2008).  See ref 1 in this paper for original reference to Qmin in
Jonsson, Mills, Jacobsen.

.. _Bitzek:

**(Bitzek)** Bitzek, Koskinen, Gahler, Moseler, Gumbsch, Phys Rev Lett,
97, 170201 (2006).

.. _Guenole:

**(Guenole)** Guenole, Noehring, Vaid, Houlle, Xie, Prakash, Bitzek,
Comput Mater Sci, (2020), in press (arXiv:190802038).
