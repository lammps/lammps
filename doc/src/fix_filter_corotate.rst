.. index:: fix filter/corotate

fix filter/corotate command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID filter/corotate keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* one or more constraint/value pairs are appended
* constraint = *b* or *a* or *t* or *m*

  .. parsed-literal::

       *b* values = one or more bond types
       *a* values = one or more angle types
       *t* values = one or more atom types
       *m* value = one or more mass values

Examples
""""""""

.. code-block:: LAMMPS

   timestep 8
   run_style respa 3 2 8 bond 1 pair 2 kspace 3
   fix cor all filter/corotate m 1.0

   fix cor all filter/corotate b 4 19 a 3 5 2

Description
"""""""""""

This fix implements a corotational filter for a mollified impulse
method. In biomolecular simulations, it allows the usage of larger
timesteps for long-range electrostatic interactions.  For details, see
:ref:`(Fath) <Fath2017>`.

When using :doc:`run_style respa <run_style>` for a biomolecular
simulation with high-frequency covalent bonds, the outer time-step is
restricted to below ~ 4fs due to resonance problems. This fix filters
the outer stage of the respa and thus a larger (outer) time-step can
be used. Since in large biomolecular simulations the computation of
the long-range electrostatic contributions poses a major bottleneck,
this can significantly accelerate the simulation.

The filter computes a cluster decomposition of the molecular structure
following the criteria indicated by the options a, b, t and m. This
process is similar to the approach in :doc:`fix shake <fix_shake>`,
however, the clusters are not kept constrained. Instead, the position
is slightly modified only for the computation of long-range forces. A
good cluster decomposition constitutes in building clusters which
contain the fastest covalent bonds inside clusters.

If the clusters are chosen suitably, the :doc:`run_style respa <run_style>` is stable for outer time-steps of at least 8fs.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about these fixes is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to these fixes.  No global or per-atom quantities are
stored by these fixes for access by various :doc:`output commands <Howto_output>`.  No parameter of these fixes can be used
with the *start/stop* keywords of the :doc:`run <run>` command.  These
fixes are not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the USER-MISC package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

Currently, it does not support :doc:`molecule templates <molecule>`.

Related commands
""""""""""""""""

**Default:** none

----------

.. _Fath2017:

**(Fath)** Fath, Hochbruck, Singh, J Comp Phys, 333, 180-198 (2017).
