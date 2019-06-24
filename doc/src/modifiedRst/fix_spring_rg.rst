.. index:: fix spring/rg

fix spring/rg command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID spring/rg K RG0

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring/rg = style name of this fix command
* K = harmonic force constant (force/distance units)
* RG0 = target radius of gyration to constrain to (distance units)


.. parsed-literal::

     if RG0 = NULL, use the current RG as the target value

Examples
""""""""


.. parsed-literal::

   fix 1 protein spring/rg 5.0 10.0
   fix 2 micelle spring/rg 5.0 NULL

Description
"""""""""""

Apply a harmonic restraining force to atoms in the group to affect
their central moment about the center of mass (radius of gyration).
This fix is useful to encourage a protein or polymer to fold/unfold
and also when sampling along the radius of gyration as a reaction
coordinate (i.e. for protein folding).

The radius of gyration is defined as RG in the first formula.  The
energy of the constraint and associated force on each atom is given by
the second and third formulas, when the group is at a different RG
than the target value RG0.

.. math source doc: src/Eqs/fix_spring_rg.tex
.. math::

   {R_G}^2 = \frac{1}{M}\sum_{i}^{N}{m_{i}\left( x_{i} -
   \frac{1}{M}\sum_{j}^{N}{m_{j}x_{j}} \right)^{2}}

>>>image was here
   E = K\left( R_G - R_{G0} \right)^{2}

>>>image was here
   F_{i} = 2K\frac{m_{i}}{M}\left( 1-\frac{R_{G0}}{R_G}
   \right)\left( x_{i} - \frac{1}{M}\sum_{j}^{N}{m_{j}x_{j}} \right)


The (xi - center-of-mass) term is computed taking into account
periodic boundary conditions, m\_i is the mass of the atom, and M is
the mass of the entire group.  Note that K is thus a force constant
for the aggregate force on the group of atoms, not a per-atom force.

If RG0 is specified as NULL, then the RG of the group is computed at
the time the fix is specified, and that value is used as the target.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

The :doc:`fix\_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix spring <fix_spring>`, :doc:`fix spring/self <fix_spring_self>`
:doc:`fix drag <fix_drag>`, :doc:`fix smd <fix_smd>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
