.. index:: fix nvt/manifold/rattle

fix nvt/manifold/rattle command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nvt/manifold/rattle tol maxit manifold manifold-args keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvt/manifold/rattle = style name of this fix command
* tol = tolerance to which Newton iteration must converge
* maxit = maximum number of iterations to perform
* manifold = name of the manifold
* manifold-args = parameters for the manifold
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *temp* or *tchain* or *every*
       *temp* values = Tstart Tstop Tdamp
         Tstart, Tstop = external temperature at start/end of run
         Tdamp = temperature damping parameter (time units)
       *tchain* value = N
         N = length of thermostat chain (1 = single thermostat)
       *every* value = N
         N = print info about iteration every N steps. N = 0 means no output

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/manifold/rattle 1e-4 10 cylinder 3.0 temp 1.0 1.0 10.0

Description
"""""""""""

This fix combines the RATTLE-based :ref:`(Andersen) <Andersen2>` time
integrator of :doc:`fix nve/manifold/rattle <fix_nve_manifold_rattle>`
:ref:`(Paquay) <Paquay3>` with a Nose-Hoover-chain thermostat to sample the
canonical ensemble of particles constrained to a curved surface
(manifold). This sampling does suffer from discretization bias of
O(dt).  For a list of currently supported manifolds and their
parameters, see the :doc:`Howto manifold <Howto_manifold>` doc page.

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
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

----------

Related commands
""""""""""""""""

:doc:`fix nve/manifold/rattle <fix_nvt_manifold_rattle>`, :doc:`fix manifoldforce <fix_manifoldforce>` **Default:** every = 0

----------

.. _Andersen2:

**(Andersen)** Andersen, J. Comp. Phys. 52, 24, (1983).

.. _Paquay3:

**(Paquay)** Paquay and Kusters, Biophys. J., 110, 6, (2016).
preprint available at `arXiv:1411.3019 <https://arxiv.org/abs/1411.3019/>`_.
