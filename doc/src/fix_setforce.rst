.. index:: fix setforce

fix setforce command
====================

fix setforce/kk command
=======================

fix setforce/spin command
=========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID setforce fx fy fz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* setforce = style name of this fix command
* fx,fy,fz = force component values
* any of fx,fy,fz can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to have added force

Examples
""""""""

.. parsed-literal::

   fix freeze indenter setforce 0.0 0.0 0.0
   fix 2 edge setforce NULL 0.0 0.0
   fix 1 edge setforce/spin 0.0 0.0 0.0
   fix 2 edge setforce NULL 0.0 v_oscillate

Description
"""""""""""

Set each component of force on each atom in the group to the specified
values fx,fy,fz.  This erases all previously computed forces on the
atom, though additional fixes could add new forces.  This command can
be used to freeze certain atoms in the simulation by zeroing their
force, either for running dynamics or performing an energy
minimization.  For dynamics, this assumes their initial velocity is
also zero.

Any of the fx,fy,fz values can be specified as NULL which means do not
alter the force component in that dimension.

Any of the 3 quantities defining the force components can be specified
as an equal-style or atom-style :doc:`variable <variable>`, namely *fx*\ ,
*fy*\ , *fz*\ .  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
force component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent force field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent force
field with optional time-dependence as well.

If the *region* keyword is used, the atom must also be in the
specified geometric :doc:`region <region>` in order to have force added
to it.

----------

Style *spin* suffix sets the components of the magnetic precession
vectors instead of the mechanical forces. This also erases all
previously computed magnetic precession vectors on the atom, though
additional magnetic fixes could add new forces.

This command can be used to freeze the magnetic moment of certain
atoms in the simulation by zeroing their precession vector.

All options defined above remain valid, they just apply to the magnetic
precession vectors instead of the forces.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

The region keyword is also supported by Kokkos, but a Kokkos-enabled
region must be used. See the region :doc:`region <region>` command for
more information.

These accelerated styles are part of the r Kokkos package.  They are
only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by
this fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is setting the forces to the desired values; on all
other levels, the force is set to 0.0 for the atoms in the fix group,
so that setforce values are not counted multiple times. Default is to
to override forces at the outermost level.

This fix computes a global 3-vector of forces, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
force on the group of atoms before the forces on individual atoms are
changed by the fix.  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command, but you cannot set
forces to any value besides zero when performing a minimization.  Use
the :doc:`fix addforce <fix_addforce>` command if you want to apply a
non-zero force to atoms during a minimization.

Restrictions
""""""""""""

The fix *setforce/spin* only makes sense when LAMMPS was built with the
SPIN package.

Related commands
""""""""""""""""

:doc:`fix addforce <fix_addforce>`, :doc:`fix aveforce <fix_aveforce>`

**Default:** none
