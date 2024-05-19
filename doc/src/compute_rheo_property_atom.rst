.. index:: compute rheo/property/atom

compute rheo/property/atom command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID rheo/property/atom input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* rheo/property/atom = style name of this compute command
* input = one or more atom attributes

  .. parsed-literal::

       possible attributes = phase, surface, surface/r,
                             surface/divr, surface/n/x, surface/n/y,
                             surface/n/z, coordination, cv, shift/v/x,
                             shift/v/y, shift/v/z, energy, temperature, heatflow,
                             conductivity, cv, viscosity, pressure,
                             rho, grad/v/xx, grad/v/xy, grad/v/xz,
                             grad/v/yx, grad/v/yy/, grad/v/yz, grad/v/zx,
                             grad/v/zy, grad/v/zz

  .. parsed-literal::

           *phase* = atom phase state
           *surface* = atom surface status
           *surface/r* = atom distance from the surface
           *surface/divr* = divergence of position at atom position
           *surface/n/\** = surface normal vector
           *coordination* = coordination number
           *shift/v/\** = atom shifting velocity
           *energy* = atom energy
           *temperature* = atom temperature
           *heatflow* = atom heat flow
           *conductivity* = atom conductivity
           *cv* = atom specific heat
           *viscosity* = atom viscosity
           *pressure* = atom pressure
           *rho* = atom density
           *grad/v/\** = atom velocity gradient
           *nbond/shell* = number of oxide bonds

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all rheo/property/atom phase surface/r pressure
   compute 2 all rheo/property/atom shift/v/x grad/v/xx

Description
"""""""""""

.. versionadded:: TBD

Define a computation that stores atom attributes specific to the RHEO
package for each atom in the group.  This is useful so that the values
can be used by other :doc:`output commands <Howto_output>` that take
computes as inputs. See for example, the
:doc:`compute reduce <compute_reduce>`,
:doc:`fix ave/atom <fix_ave_atom>`,
:doc:`fix ave/histo <fix_ave_histo>`,
:doc:`fix ave/chunk <fix_ave_chunk>`, and
:doc:`atom-style variable <variable>` commands.

Many properties require their respective fixes, listed below in related
commands, be defined. For instance, the *viscosity* attribute is the
viscosity of a particle calculated by
:doc:`fix rheo/viscous <fix_rheo_viscosity>`. The meaning of less obvious
properties is described below.

The *phase* property indicates whether the particle is in a fluid state,
a value of 0, or a solid state, a value of 1.

The *surface* property indicates the surface designation produced by
the *interface/reconstruct* option of :doc:`fix rheo <fix_rheo>`. Bulk
particles have a value of 0, surface particles have a value of 1, and
splash particles have a value of 2. The *surface/r* property is the
distance from the surface, up to the kernel cutoff length. Surface particles
have a value of 0. The *surface/n* properties are the components of the
surface normal vector.

The *shift/v* properties are the components of the shifting velocity
produced by the *shift* option of :doc:`fix rheo <fix_rheo>`.

The *surface/n/\** and *shift/v/\** attributes are vectors that require
specification of the *x*, *y*, or *z* component, e.g. *surface/n/x*.

The *nbond/shell* property is the number of shell bonds that have been
activated from :doc:`bond style rheo/shell <bond_rheo_shell>`.

The *grad/v/\** attribute is a tensor and requires specification of
the *xx*, *yy*, *zz*, *xy*, *xz*, *yx*, *yz*, *zx*, or *zy* component,
e.g. *grad/v/xy*.

The values are stored in a per-atom vector or array as discussed
below.  Zeroes are stored for atoms not in the specified group or for
quantities that are not defined for a particular particle in the group

Output info
"""""""""""

This compute calculates a per-atom vector or per-atom array depending
on the number of input values.  If a single input is specified, a
per-atom vector is produced.  If two or more inputs are specified, a
per-atom array is produced where the number of columns = the number of
inputs.  The vector or array can be accessed by any command that uses
per-atom values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

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
