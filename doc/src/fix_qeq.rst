.. index:: fix qeq/point

fix qeq/point command
=====================

fix qeq/shielded command
========================

fix qeq/slater command
======================

fix qeq/dynamic command
=======================

fix qeq/fire command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID style Nevery cutoff tolerance maxiter qfile keyword ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *qeq/point* or *qeq/shielded* or *qeq/slater* or *qeq/dynamic* or *qeq/fire*
* Nevery = perform charge equilibration every this many steps
* cutoff = global cutoff for charge-charge interactions (distance unit)
* tolerance = precision to which charges will be equilibrated
* maxiter = maximum iterations to perform charge equilibration
* qfile = a filename with QEq parameters or *coul/streitz* or *reax/c*
* zero or more keyword/value pairs may be appended
* keyword = *alpha* or *qdamp* or *qstep*

  .. parsed-literal::

       *alpha* value = Slater type orbital exponent (qeq/slater only)
       *qdamp* value = damping factor for damped dynamics charge solver (qeq/dynamic and qeq/fire only)
       *qstep* value = time step size for damped dynamics charge solver (qeq/dynamic and qeq/fire only)



Examples
""""""""


.. parsed-literal::

   fix 1 all qeq/point 1 10 1.0e-6 200 param.qeq1
   fix 1 qeq qeq/shielded 1 8 1.0e-6 100 param.qeq2
   fix 1 all qeq/slater 5 10 1.0e-6 100 params alpha 0.2
   fix 1 qeq qeq/dynamic 1 12 1.0e-3 100 my_qeq
   fix 1 all qeq/fire 1 10 1.0e-3 100 my_qeq qdamp 0.2 qstep 0.1

Description
"""""""""""

Perform the charge equilibration (QEq) method as described in :ref:`(Rappe and Goddard) <Rappe1>` and formulated in :ref:`(Nakano) <Nakano1>` (also known
as the matrix inversion method) and in :ref:`(Rick and Stuart) <Rick1>` (also
known as the extended Lagrangian method) based on the
electronegativity equilization principle.

These fixes can be used with any :doc:`pair style <pair_style>` in
LAMMPS, so long as per-atom charges are defined.  The most typical
use-case is in conjunction with a :doc:`pair style <pair_style>` that
performs charge equilibration periodically (e.g. every timestep), such
as the ReaxFF or Streitz-Mintmire potential.
But these fixes can also be used with
potentials that normally assume per-atom charges are fixed, e.g. a
:doc:`Buckingham <pair_buck>` or :doc:`LJ/Coulombic <pair_lj>` potential.

Because the charge equilibration calculation is effectively
independent of the pair style, these fixes can also be used to perform
a one-time assignment of charges to atoms.  For example, you could
define the QEq fix, perform a zero-timestep run via the :doc:`run <run>`
command without any pair style defined which would set per-atom
charges (based on the current atom configuration), then remove the fix
via the :doc:`unfix <unfix>` command before performing further dynamics.

.. note::

   Computing and using charge values different from published
   values defined for a fixed-charge potential like Buckingham or CHARMM
   or AMBER, can have a strong effect on energies and forces, and
   produces a different model than the published versions.

.. note::

   The :doc:`fix qeq/comb <fix_qeq_comb>` command must still be used
   to perform charge equilibration with the :doc:`COMB potential <pair_comb>`.  The :doc:`fix qeq/reax <fix_qeq_reax>`
   command can be used to perform charge equilibration with the :doc:`ReaxFF force field <pair_reaxc>`, although fix qeq/shielded yields the
   same results as fix qeq/reax if *Nevery*\ , *cutoff*\ , and *tolerance*
   are the same.  Eventually the fix qeq/reax command will be deprecated.

The QEq method minimizes the electrostatic energy of the system (or
equalizes the derivative of energy with respect to charge of all the
atoms) by adjusting the partial charge on individual atoms based on
interactions with their neighbors within *cutoff*\ .  It requires a few
parameters, in *metal* units, for each atom type which provided in a
file specified by *qfile*\ .  The file has the following format


.. parsed-literal::

   1 chi eta gamma zeta qcore
   2 chi eta gamma zeta qcore
   ...
   Ntype chi eta gamma zeta qcore

There have to be parameters given for every atom type. Wildcard entries
are possible using the same syntax as elsewhere in LAMMPS
(i.e., n\*m, n\*, \*m, \*). Later entries will overwrite previous ones.
Empty lines or any text following the pound sign (#) are ignored.
Each line starts with the atom type followed by five parameters.
Only a subset of the parameters is used by each QEq style as described
below, thus the others can be set to 0.0 if desired, but all five
entries per line are required.

* *chi* = electronegativity in energy units
* *eta* = self-Coulomb potential in energy units
* *gamma* = shielded Coulomb constant defined by :ref:`ReaxFF force field <vanDuin>` in distance units
* *zeta* = Slater type orbital exponent defined by the :ref:`Streitz-Mintmire <Streitz1>` potential in reverse distance units
* *qcore* = charge of the nucleus defined by the :ref:`Streitz-Mintmire potential <Streitz1>` potential in charge units

The *qeq/point* style describes partial charges on atoms as point
charges.  Interaction between a pair of charged particles is 1/r,
which is the simplest description of the interaction between charges.
Only the *chi* and *eta* parameters from the *qfile* file are used.
Note that Coulomb catastrophe can occur if repulsion between the pair
of charged particles is too weak.  This style solves partial charges
on atoms via the matrix inversion method.  A tolerance of 1.0e-6 is
usually a good number.

The *qeq/shielded* style describes partial charges on atoms also as
point charges, but uses a shielded Coulomb potential to describe the
interaction between a pair of charged particles.  Interaction through
the shielded Coulomb is given by equation (13) of the :ref:`ReaxFF force field <vanDuin>` paper.  The shielding accounts for charge overlap
between charged particles at small separation.  This style is the same
as :doc:`fix qeq/reax <fix_qeq_reax>`, and can be used with :doc:`pair_style reax/c <pair_reaxc>`.  Only the *chi*\ , *eta*\ , and *gamma*
parameters from the *qfile* file are used. When using the string
*reax/c* as filename, these parameters are extracted directly from
an active *reax/c* pair style.  This style solves partial
charges on atoms via the matrix inversion method.  A tolerance of
1.0e-6 is usually a good number.

The *qeq/slater* style describes partial charges on atoms as spherical
charge densities centered around atoms via the Slater 1\ *s* orbital, so
that the interaction between a pair of charged particles is the
product of two Slater 1\ *s* orbitals.  The expression for the Slater
1\ *s* orbital is given under equation (6) of the
:ref:`Streitz-Mintmire <Streitz1>` paper.  Only the *chi*\ , *eta*\ , *zeta*\ , and
*qcore* parameters from the *qfile* file are used. When using the string
*coul/streitz* as filename, these parameters are extracted directly from
an active *coul/streitz* pair style.  This style solves
partial charges on atoms via the matrix inversion method.  A tolerance
of 1.0e-6 is usually a good number.  Keyword *alpha* can be used to
change the Slater type orbital exponent.

The *qeq/dynamic* style describes partial charges on atoms as point
charges that interact through 1/r, but the extended Lagrangian method
is used to solve partial charges on atoms.  Only the *chi* and *eta*
parameters from the *qfile* file are used.  Note that Coulomb
catastrophe can occur if repulsion between the pair of charged
particles is too weak.  A tolerance of 1.0e-3 is usually a good
number.  Keyword *qdamp* can be used to change the damping factor, while
keyword *qstep* can be used to change the time step size.

The :ref:`\ *qeq/fire*\ <Shan>` style describes the same charge model and charge
solver as the *qeq/dynamic* style, but employs a FIRE minimization
algorithm to solve for equilibrium charges.
Keyword *qdamp* can be used to change the damping factor, while
keyword *qstep* can be used to change the time step size.

Note that *qeq/point*\ , *qeq/shielded*\ , and *qeq/slater* describe
different charge models, whereas the matrix inversion method and the
extended Lagrangian method (\ *qeq/dynamic* and *qeq/fire*\ ) are
different solvers.

Note that *qeq/point*\ , *qeq/dynamic* and *qeq/fire* styles all describe
charges as point charges that interact through 1/r relationship, but
solve partial charges on atoms using different solvers.  These three
styles should yield comparable results if
the QEq parameters and *Nevery*\ , *cutoff*\ , and *tolerance* are the
same.  Style *qeq/point* is typically faster, *qeq/dynamic* scales
better on larger sizes, and *qeq/fire* is faster than *qeq/dynamic*\ .

.. note::

   To avoid the evaluation of the derivative of charge with respect
   to position, which is typically ill-defined, the system should have a
   zero net charge.

.. note::

   Developing QEq parameters (chi, eta, gamma, zeta, and qcore) is
   non-trivial.  Charges on atoms are not guaranteed to equilibrate with
   arbitrary choices of these parameters.  We do not develop these QEq
   parameters.  See the examples/qeq directory for some examples.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about these fixes is written to :doc:`binary restart files <restart>`.  No global scalar or vector or per-atom
quantities are stored by these fixes for access by various :doc:`output commands <Howto_output>`.  No parameter of these fixes can be used
with the *start/stop* keywords of the :doc:`run <run>` command.

Thexe fixes are invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


These fixes are part of the QEQ package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix qeq/reax <fix_qeq_reax>`, :doc:`fix qeq/comb <fix_qeq_comb>`

**Default:** none


----------


.. _Rappe1:



**(Rappe and Goddard)** A. K. Rappe and W. A. Goddard III, J Physical
Chemistry, 95, 3358-3363 (1991).

.. _Nakano1:



**(Nakano)** A. Nakano, Computer Physics Communications, 104, 59-69 (1997).

.. _Rick1:



**(Rick and Stuart)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chemical Physics
101, 16141 (1994).

.. _Streitz1:



**(Streitz-Mintmire)** F. H. Streitz, J. W. Mintmire, Physical Review B, 50,
16, 11996 (1994)

.. _vanDuin:



**(ReaxFF)** A. C. T. van Duin, S. Dasgupta, F. Lorant, W. A. Goddard III, J
Physical Chemistry, 105, 9396-9049 (2001)

.. _Shan:



**(QEq/Fire)** T.-R. Shan, A. P. Thompson, S. J. Plimpton, in preparation
