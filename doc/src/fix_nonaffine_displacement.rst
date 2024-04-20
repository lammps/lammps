.. index:: fix nonaffine/displacement

fix nonaffine/displacement command
==================================

Syntax
""""""

.. parsed-literal::

   fix ID group nonaffine/displacement style args reference/style nstep

* ID, group are documented in :doc:`fix <fix>` command
* nonaffine/displacement = style name of this fix command
* nevery = calculate nonaffine displacement every this many timesteps
* style = *d2min* or *integrated*

  .. parsed-literal::

       *d2min* args = cutoff args
         cutoff = *type* or *radius* or *custom*
           *type* args = none, cutoffs determined by atom types
           *radius* args = none, cutoffs determined based on atom diameters (atom style sphere)
           *custom* args = *rmax*, cutoff set by a constant numeric value *rmax* (distance units)
       *integrated* args = none

* reference/style = *fixed* or *update* or *offset*

  .. parsed-literal::

       *fixed* = use a fixed reference frame at *nstep*
       *update* = update the reference frame every *nstep* timesteps
       *offset* = update the reference frame *nstep* timesteps before calculating the nonaffine displacement

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nonaffine/displacement 100 integrated update 100
   fix 1 all nonaffine/displacement 1000 d2min type fixed 0
   fix 1 all nonaffine/displacement 1000 d2min custom 2.0 offset 100

Description
"""""""""""

.. versionadded:: 7Feb2024

This fix computes different metrics of the nonaffine displacement of
particles. The first metric, *d2min* calculates the :math:`D^2_\mathrm{min}`
nonaffine displacement by Falk and Langer in :ref:`(Falk) <d2min-Falk>`.
For each atom, the fix computes the two tensors

.. math::

   X = \sum_{\mathrm{neighbors}} \vec{r} \left(\vec{r}_{0} \right)^T

and

.. math::

   Y = \sum_{\mathrm{neighbors}} \vec{r}_0 \left(\vec{r}_{0} \right)^T

where the neighbors include all other atoms within the distance criterion
set by the cutoff option, discussed below, :math:`\vec{r}` is the current
displacement between particles, and :math:`\vec{r}_0` is the reference
displacement. A deformation gradient tensor is then calculated as
:math:`F = X Y^{-1}` from which

.. math::

    D^2_\mathrm{min} = \sum_{\mathrm{neighbors}} \left| \vec{r} - F \vec{r}_0 \right|^2

and a strain tensor is calculated :math:`E = F F^{T} - I` where :math:`I`
is the identity tensor. This calculation is only performed on timesteps that
are a multiple of *nevery* (including timestep zero). Data accessed before
this occurs will simply be zeroed.

The *integrated* style simply integrates the velocity of particles
every timestep to calculate a displacement. This style only works if
used in conjunction with another fix that deforms the box and displaces
atom positions such as :doc:`fix deform <fix_deform>` with remap x,
:doc:`fix press/berendsen <fix_press_berendsen>`, or :doc:`fix nh <fix_nh>`.

Both of these methods require defining a reference state. With the *fixed* reference
style, the user picks a specific timestep *nstep* at which particle positions are saved.
If peratom data is accessed from this compute prior to this timestep, it will simply be
zeroed. The *update* reference style implies the reference state will be updated every
*nstep* timesteps. The *offset* reference will update the reference state *nstep*
timesteps before a multiple of *nevery* timesteps.


----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The reference state is saved to :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a peratom array with 3 columns, which can be accessed
by indices 1-3 using any command that uses per-atom values from a fix
as input.

For the *integrated* style, the three columns are the nonaffine
displacements in the x, y, and z directions. For the *d2min* style,
the three columns are the calculated :math:`\sqrt{D^2_\mathrm{min}}`, the
volumetric strain, and the deviatoric strain.

Restrictions
""""""""""""

This compute is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

none

Default
"""""""

none

----------

.. _d2min-Falk:

**(Falk)** Falk and Langer PRE, 57, 7192 (1998).
