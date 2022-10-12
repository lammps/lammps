.. index:: fix nve/manifold/rattle

fix nve/manifold/rattle command
===============================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nve/manifold/rattle tol maxit manifold manifold-args keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/manifold/rattle = style name of this fix command
* tol = tolerance to which Newton iteration must converge
* maxit = maximum number of iterations to perform
* manifold = name of the manifold
* manifold-args = parameters for the manifold
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *every*
       *every* values = N
         N = print info about iteration every N steps. N = 0 means no output

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/manifold/rattle 1e-4 10 sphere 5.0
   fix step all nve/manifold/rattle 1e-8 100 ellipsoid 2.5 2.5 5.0 every 25

Description
"""""""""""

Perform constant NVE integration to update position and velocity for
atoms constrained to a curved surface (manifold) in the group each
timestep. The constraint is handled by RATTLE :ref:`(Andersen) <Andersen1>`
written out for the special case of single-particle constraints as
explained in :ref:`(Paquay) <Paquay2>`.  V is volume; E is energy. This way,
the dynamics of particles constrained to curved surfaces can be
studied. If combined with :doc:`fix langevin <fix_langevin>`, this
generates Brownian motion of particles constrained to a curved
surface. For a list of currently supported manifolds and their
parameters, see the :doc:`Howto manifold <Howto_manifold>` doc page.

Note that the particles must initially be close to the manifold in
question. If not, RATTLE will not be able to iterate until the
constraint is satisfied, and an error is generated. For simple
manifolds this can be achieved with *region* and *create_atoms*
commands, but for more complex surfaces it might be more useful to
write a script.

The manifold args may be equal-style variables, like so:

.. code-block:: LAMMPS

   variable R equal "ramp(5.0,3.0)"
   fix shrink_sphere all nve/manifold/rattle 1e-4 10 sphere v_R

In this case, the manifold parameter will change in time according to
the variable.  This is not a problem for the time integrator as long
as the change of the manifold is slow with respect to the dynamics of
the particles.  Note that if the manifold has to exert work on the
particles because of these changes, the total energy might not be
conserved.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix is part of the MANIFOLD package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

----------

Related commands
""""""""""""""""

:doc:`fix nvt/manifold/rattle <fix_nvt_manifold_rattle>`, :doc:`fix manifoldforce <fix_manifoldforce>`

Default
"""""""

every = 0, tchain = 3


----------

.. _Andersen1:

**(Andersen)** Andersen, J. Comp. Phys. 52, 24, (1983).

.. _Paquay2:

**(Paquay)** Paquay and Kusters, Biophys. J., 110, 6, (2016).
preprint available at `arXiv:1411.3019 <https://arxiv.org/abs/1411.3019/>`_.
