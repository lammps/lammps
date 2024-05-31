.. index:: fix deform/pressure

fix deform/pressure command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID deform/pressure N parameter style args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* deform/pressure = style name of this fix command
* N = perform box deformation every this many timesteps
* one or more parameter/arg sequences may be appended

  .. parsed-literal::

     parameter = *x* or *y* or *z* or *xy* or *xz* or *yz* or *box*
       *x*, *y*, *z* args = style value(s)
         style = *final* or *delta* or *scale* or *vel* or *erate* or *trate* or *volume* or *wiggle* or *variable* or *pressure* or *pressure/mean*
           *pressure* values = target gain
             target = target pressure (pressure units)
             gain = proportional gain constant (1/(time * pressure) or 1/time units)
           *pressure/mean* values = target gain
             target = target pressure (pressure units)
             gain = proportional gain constant (1/(time * pressure) or 1/time units)
           NOTE: All other styles are documented by the :doc:`fix deform <fix_deform>` command

       *xy*, *xz*, *yz* args = style value
         style = *final* or *delta* or *vel* or *erate* or *trate* or *wiggle* or *variable* or *pressure* or *erate/rescale*
           *pressure* values = target gain
             target = target pressure (pressure units)
             gain = proportional gain constant (1/(time * pressure) or 1/time units)
           *erate/rescale* value = R
             R = engineering shear strain rate (1/time units)
           NOTE: All other styles are documented by the :doc:`fix deform <fix_deform>` command

       *box* = style value
         style = *volume* or *pressure*
           *volume* value = none = isotropically adjust system to preserve volume of system
           *pressure* values = target gain
             target = target mean pressure (pressure units)
             gain = proportional gain constant (1/(time * pressure) or 1/time units)

* zero or more keyword/value pairs may be appended
* keyword = *remap* or *flip* or *units* or *couple* or *vol/balance/p* or *max/rate* or *normalize/pressure*

  .. parsed-literal::

       *couple* value = *none* or *xyz* or *xy* or *yz* or *xz*
         couple pressure values of various dimensions
       *vol/balance/p* value = *yes* or *no*
         Modifies the behavior of the *volume* option to try and balance pressures
       *max/rate* value = *rate*
         rate = maximum strain rate for pressure control
       *normalize/pressure* value = *yes* or *no*
         Modifies pressure controls such that the deviation in pressure is normalized by the target pressure
       NOTE: All other keywords are documented by the :doc:`fix deform <fix_deform>` command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all deform/pressure 1 x pressure 2.0 0.1 normalize/pressure yes max/rate 0.001
   fix 1 all deform/pressure 1 x trate 0.1 y volume z volume vol/balance/p yes
   fix 1 all deform/pressure 1 x trate 0.1 y pressure/mean 0.0 1.0 z pressure/mean 0.0 1.0

Description
"""""""""""

.. versionadded:: 17Apr2024

This fix is an extension of the :doc:`fix deform <fix_deform>`
command, which allows all of its options to be used as well as new
pressure-based controls implemented by this command.

All arguments described on the :doc:`fix deform <fix_deform>` doc page
also apply to this fix unless otherwise noted below.  The rest of this
doc page explains the arguments specific to this fix.  Note that a
simulation can define only a single deformation command: fix deform or
fix deform/pressure.

----------

For the *x*, *y*, and *z* parameters, this is the meaning of the
styles and values provided by this fix.

The *pressure* style adjusts a dimension's box length to control the
corresponding component of the pressure tensor. This option attempts to
maintain a specified target pressure using a linear controller where the
box length :math:`L` evolves according to the equation

.. math::

   \frac{d L(t)}{dt} = L(t) k (P_t - P)

where :math:`k` is a proportional gain constant, :math:`P_t` is the target
pressure, and :math:`P` is the current pressure along that dimension. This
approach is similar to the method used to control the pressure by
:doc:`fix press/berendsen <fix_press_berendsen>`. The target pressure
accepts either a constant numeric value or a LAMMPS :ref:`variable <variable>`.
Notably, this variable can be a function of time or other components of
the pressure tensor. By default, :math:`k` has units of 1/(time * pressure)
although this will change if the *normalize/pressure* option is set as
:ref:`discussed below <deform_normalize>`. There is no proven method
to choosing an appropriate value of :math:`k` as it will depend on the
specific details of a simulation. Testing different values is recommended.

By default, there is no limit on the resulting strain rate in any dimension.
A maximum limit can be applied using the :ref:`max/rate <deform_max_rate>`
option. Akin to :doc:`fix nh <fix_nh>`, pressures in different dimensions
can be coupled using the :ref:`couple <deform_couple>` option. This means
the instantaneous pressure along coupled dimensions are averaged and the box
strains identically along the coupled dimensions.

The *pressure/mean* style changes a dimension's box length to maintain
a constant mean pressure defined as the trace of the pressure tensor.
This option has identical arguments to the *pressure* style and a similar
functional equation, except the current and target pressures refer to the
mean trace of the pressure tensor. All options for the *pressure* style
also apply to the *pressure/mean* style except for the
:ref:`couple <deform_couple>` option.

Note that while this style can be identical to coupled *pressure* styles,
it is generally not the same. For instance in 2D, a coupled *pressure*
style in the *x* and *y* dimensions would be equivalent to using the
*pressure/mean* style with identical settings in each dimension. However,
it would not be the same if settings (e.g. gain constants) were used in
the *x* and *y* dimensions or if the *pressure/mean* command was only applied
along one dimension.

----------

For the *xy*, *xz*, and *yz* parameters, this is the meaning of the
styles and values provided by this fix.  Note that changing the
tilt factors of a triclinic box does not change its volume.

The *pressure* style adjusts a tilt factor to control the corresponding
off-diagonal component of the pressure tensor. This option attempts to
maintain a specified target value using a linear controller where the
tilt factor T evolves according to the equation

.. parsed-literal::

   \frac{d T(t)}{dt} = L(t) k (P - P_t)

where :math:`k` is a proportional gain constant, :math:`P_t` is the
target pressure, :math:`P` is the current pressure, and :math:`L` is
the perpendicular box length. The target pressure accepts either a
constant numeric value or a LAMMPS :ref:`variable
<variable>`. Notably, this variable can be a function of time or other
components of the pressure tensor. By default, :math:`k` has units of
1/(time * pressure) although this will change if the
*normalize/pessure* option is set as :ref:`discussed below
<deform_normalize>`.  There is no proven method to choosing an
appropriate value of :math:`k` as it will depend on the specific
details of a simulation and testing different values is
recommended. One can also apply a maximum limit to the magnitude of
the applied strain using the :ref:`max/rate <deform_max_rate>` option.

The *erate/rescale* style operates similarly to the *erate* style with
a specified strain rate in units of 1/time. The difference is that
the change in the tilt factor will depend on the current length of
the box perpendicular to the shear direction, L, instead of the
original length, L0. The tilt factor T as a function of time will
change as

.. parsed-literal::

   T(t) = T(t-1) + L\*erate\* \Delta t

where T(t-1) is the tilt factor on the previous timestep and :math:`\Delta t`
is the timestep size. This option may be useful in scenarios where
L changes in time.

----------

The *box* parameter provides an additional control over the *x*, *y*,
and *z* box lengths by isotropically dilating or contracting the box
to either maintain a fixed mean pressure or volume. This isotropic
scaling is applied after the box is deformed by the above *x*, *y*,
*z*, *xy*, *xz*, and *yz* styles, acting as a second deformation
step. This parameter will change the overall strain rate in the *x*,
*y*, or *z* dimensions.  This parameter can only be used in
combination with the *x*, *y*, or *z* commands: *vel*, *erate*,
*trate*, *pressure*, or *wiggle*. This is the meaning of its styles
and values.

The *volume* style isotropically scales box lengths to maintain a constant
box volume in response to deformation from other parameters. This style
may be useful in scenarios where one wants to apply a constant deviatoric
pressure using *pressure* styles in the *x*, *y*, and *z* dimensions (
deforming the shape of the box), while maintaining a constant volume.

The *pressure* style isotropically scales box lengths in an attempt to
maintain a target mean pressure (the trace of the pressure tensor) of the
system. This is accomplished by isotropically scaling all box lengths
:math:`L` by an additional factor of :math:`k (P_t - P_m)` where :math:`k`
is the proportional gain constant, :math:`P_t` is the target pressure, and
:math:`P_m` is the current mean pressure. This style may be useful in
scenarios where one wants to apply a constant deviatoric strain rate
using various strain-based styles (e.g. *trate*) along the *x*, *y*, and *z*
dimensions (deforming the shape of the box), while maintaining a mean pressure.

----------

The optional keywords provided by this fix are described below.

.. _deform_normalize:

The *normalize/pressure* keyword changes how box dimensions evolve when
using the *pressure* or *pressure/mean* deformation styles. If the
*deform/normalize* value is set to *yes*, then the deviation from the
target pressure is normalized by the absolute value of the target
pressure such that the proportional gain constant scales a percentage
error and has units of 1/time. If the target pressure is ever zero, this
will produce an error unless the *max/rate* keyword is defined,
described below, which will cap the divergence.

.. _deform_max_rate:

The *max/rate* keyword sets an upper threshold, *rate*, that limits the
maximum magnitude of the instantaneous strain rate applied in any dimension.
This keyword only applies to the *pressure* and *pressure/mean* options. If
a pressure-controlled rate is used for both *box* and either *x*, *y*, or
*z*, then this threshold will apply separately to each individual controller
such that the cumulative strain rate on a box dimension may be up to twice
the value of *rate*.

.. _deform_couple:

The *couple* keyword allows two or three of the diagonal components of
the pressure tensor to be "coupled" together for the *pressure* option.
The value specified with the keyword determines which are coupled. For
example, *xz* means the *Pxx* and *Pzz* components of the stress tensor
are coupled. *Xyz* means all 3 diagonal components are coupled. Coupling
means two things: the instantaneous stress will be computed as an average
of the corresponding diagonal components, and the coupled box dimensions
will be changed together in lockstep, meaning coupled dimensions will be
dilated or contracted by the same percentage every timestep. If a *pressure*
style is defined for more than one coupled dimension, the target pressures
and gain constants must be identical. Alternatively, if a *pressure*
style is only defined for one of the coupled dimensions, its settings are
copied to other dimensions with undefined styles. *Couple xyz* can be used
for a 2d simulation; the *z* dimension is simply ignored.

.. _deform_balance:

The *vol/balance/p* keyword modifies the behavior of the *volume* style when
applied to two of the *x*, *y*, and *z* dimensions. Instead of straining
the two dimensions in lockstep, the two dimensions are allowed to
separately dilate or contract in a manner to maintain a constant
volume while simultaneously trying to keep the pressure along each
dimension equal using a method described in :ref:`(Huang2014) <Huang2014>`.

----------

If any pressure controls are used, this fix computes a temperature and
pressure each timestep. To do this, the fix creates its own computes
of style "temp" and "pressure", as if these commands had been issued:

.. code-block:: LAMMPS

   compute fix-ID_temp group-ID temp
   compute fix-ID_press group-ID pressure fix-ID_temp

See the :doc:`compute temp <compute_temp>` and :doc:`compute pressure
<compute_pressure>` commands for details.  Note that the IDs of the
new computes are the fix-ID + underscore + "temp" or fix_ID
+ underscore + "press", and the group for the new computes is the same
as the fix group.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID =
*thermo_temp* and *thermo_press*.  This means you can change the
attributes of this fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command or print this
temperature or pressure during thermodynamic output via the
:doc:`thermo_style custom <thermo_style>` command using the
appropriate compute-ID. It also means that changing attributes of
*thermo_temp* or *thermo_press* will have no effect on this fix.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix will restore the initial box settings from :doc:`binary
restart files <restart>`, which allows the fix to be properly continue
deformation, when using the start/stop options of the :doc:`run <run>`
command.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

If any pressure controls are used, the :doc:`fix_modify <fix_modify>`
*temp* and *press* options are supported by this fix, unlike in
:doc:`fix deform <fix_deform>`.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be
used in its temperature and pressure calculations.  If you do this,
note that the kinetic energy derived from the compute temperature
should be consistent with the virial term computed using all atoms for
the pressure.  LAMMPS will warn you if you choose to compute
temperature on a subset of atoms.

This fix can perform deformation over multiple runs, using the *start*
and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

You cannot apply x, y, or z deformations to a dimension that is
shrink-wrapped via the :doc:`boundary <boundary>` command.

You cannot apply xy, yz, or xz deformations to a second dimension (y
in xy) that is shrink-wrapped via the :doc:`boundary <boundary>`
command.

Related commands
""""""""""""""""

:doc:`fix deform <fix_deform>`, :doc:`change_box <change_box>`

Default
"""""""

The option defaults are normalize/pressure = no.

----------

.. _Huang2014:

**(Huang2014)** X. Huang, "Exploring critical-state behavior using DEM",
Doctoral dissertation, Imperial College. (2014). https://doi.org/10.25560/25316
