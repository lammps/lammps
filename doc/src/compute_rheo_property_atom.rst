.. index:: compute rheo/property/atom

compute rheo/property/atom command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID rheo/property/atom input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* rheo/property/atom = style name of this compute command
* input = one or more atom attributes

  .. parsed-literal::

       possible attributes = phase, surface, surface/r,
                             surface/divr, surface/n/a, coordination,
                             shift/v/a, energy, temperature, heatflow,
                             conductivity, cv, viscosity, pressure, rho,
                             grad/v/ab, stress/v/ab, stress/t/ab, nbond/shell

  .. parsed-literal::

           *phase* = atom phase state
           *surface* = atom surface status
           *surface/r* = atom distance from the surface
           *surface/divr* = divergence of position at atom position
           *surface/n/a* = a-component of surface normal vector
           *coordination* = coordination number
           *shift/v/a* = a-component of atom shifting velocity
           *energy* = atom energy
           *temperature* = atom temperature
           *heatflow* = atom heat flow
           *conductivity* = atom conductivity
           *cv* = atom specific heat
           *viscosity* = atom viscosity
           *pressure* = atom pressure
           *rho* = atom density
           *grad/v/ab* = ab-component of atom velocity gradient tensor
           *stress/v/ab* = ab-component of atom viscous stress tensor
           *stress/t/ab* = ab-component of atom total stress tensor (pressure and viscous)
           *nbond/shell* = number of oxide bonds

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all rheo/property/atom phase surface/r surface/n/* pressure
   compute 2 all rheo/property/atom shift/v/x grad/v/xx stress/v/*

Description
"""""""""""

.. versionadded:: 29Aug2024

Define a computation that stores atom attributes specific to the RHEO
package for each atom in the group.  This is useful so that the values
can be used by other :doc:`output commands <Howto_output>` that take
computes as inputs. See for example, the
:doc:`compute reduce <compute_reduce>`,
:doc:`fix ave/atom <fix_ave_atom>`,
:doc:`fix ave/histo <fix_ave_histo>`,
:doc:`fix ave/chunk <fix_ave_chunk>`, and
:doc:`atom-style variable <variable>` commands.

For vector attributes, e.g. *shift/v/*:math:`\alpha`, one must specify
:math:`\alpha` as the *x*, *y*, or *z* component, e.g. *shift/v/x*.
Alternatively, a wild card \* will include all components, *x* and *y* in
2D or *x*, *y*, and *z* in 3D.

For tensor attributes, e.g. *grad/v/*:math:`\alpha \beta`, one must specify
both :math:`\alpha` and :math:`\beta` as  *x*, *y*, or *z*, e.g. *grad/v/xy*.
Alternatively, a wild card \* will include all components. In 2D, this
includes *xx*, *xy*, *yx*, and *yy*. In 3D, this includes *xx*, *xy*, *xz*,
*yx*, *yy*, *yz*, *zx*, *zy*, and *zz*.

Many properties require their respective fixes, listed below in related
commands, be defined. For instance, the *viscosity* attribute is the
viscosity of a particle calculated by
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`. The meaning of less obvious
properties is described below.

The *phase* property indicates whether the particle is in a fluid state,
a value of 0, or a solid state, a value of 1.

The *surface* property indicates the surface designation produced by
the *interface/reconstruct* option of :doc:`fix rheo <fix_rheo>`. Bulk
particles have a value of 0, surface particles have a value of 1, and
splash particles have a value of 2. The *surface/r* property is the
distance from the surface, up to the kernel cutoff length. Surface particles
have a value of 0. The *surface/n/*:math:`\alpha` properties are the
components of the surface normal vector.

The *shift/v/*:math:`\alpha` properties are the components of the shifting
velocity produced by the *shift* option of :doc:`fix rheo <fix_rheo>`.

The *nbond/shell* property is the number of shell bonds that have been
activated from :doc:`bond style rheo/shell <bond_rheo_shell>`.

The values are stored in a per-atom vector or array as discussed
below.  Zeroes are stored for atoms not in the specified group or for
quantities that are not defined for a particular particle in the group

Output info
"""""""""""

This compute calculates a per-atom vector or per-atom array depending
on the number of input values.  Generally, if a single input is specified,
a per-atom vector is produced.  If two or more inputs are specified, a
per-atom array is produced where the number of columns = the number of
inputs. However, if a wild card \* is used for a vector or tensor, then
the number of inputs is considered to be incremented by the dimension or
the dimension squared, respectively. The vector or array can be accessed
by any command that uses per-atom values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The vector or array values will be in whatever :doc:`units <units>` the
corresponding attribute is in (e.g., density units for *rho*).

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump custom <dump>`, :doc:`compute reduce <compute_reduce>`,
:doc:`fix ave/atom <fix_ave_atom>`, :doc:`fix ave/chunk <fix_ave_chunk>`,
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`fix rheo/oxdiation <fix_rheo_oxidation>`,
:doc:`fix rheo <fix_rheo>`

Default
"""""""

none
