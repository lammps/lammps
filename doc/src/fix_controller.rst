.. index:: fix controller

fix controller command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID controller Nevery alpha Kp Ki Kd pvar setpoint cvar

* ID, group-ID are documented in :doc:`fix <fix>` command
* controller = style name of this fix command
* Nevery = invoke controller every this many timesteps
* alpha = coupling constant for PID equation (see units discussion below)
* Kp = proportional gain in PID equation (unitless)
* Ki = integral gain in PID equation (unitless)
* Kd = derivative gain in PID equation (unitless)
* pvar = process variable of form c\_ID, c\_ID[I], f\_ID, f\_ID[I], or v\_name

  .. parsed-literal::

       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID
       v_name = value calculated by an equal-style variable with name

* setpoint = desired value of process variable (same units as process variable)
* cvar = name of control variable

Examples
""""""""

.. parsed-literal::

   fix 1 all controller 100 1.0 0.5 0.0 0.0 c_thermo_temp 1.5 tcontrol
   fix 1 all controller 100 0.2 0.5 0 100.0 v_pxxwall 1.01325 xwall
   fix 1 all controller 10000 0.2 0.5 0 2000 v_avpe -3.785 tcontrol

Description
"""""""""""

This fix enables control of a LAMMPS simulation using a control loop
feedback mechanism known as a proportional-integral-derivative (PID)
controller.  The basic idea is to define a "process variable" which is
a quantity that can be monitored during a running simulation.  A
desired target value is chosen for the process variable.  A "control
variable" is also defined which is an adjustable attribute of the
running simulation, which the process variable will respond to.  The
PID controller continuously adjusts the control variable based on the
difference between the process variable and the target.

Here are examples of ways in which this fix can be used.  The
examples/pid directory contains a script that implements the simple
thermostat.

+-----------------------------------------+---------------------+---------------------+
| Goal                                    | process variable    | control variable    |
+-----------------------------------------+---------------------+---------------------+
| Simple thermostat                       | instantaneous T     | thermostat target T |
+-----------------------------------------+---------------------+---------------------+
| Find melting temperature                | average PE per atom | thermostat target T |
+-----------------------------------------+---------------------+---------------------+
| Control pressure in non-periodic system | force on wall       | position of wall    |
+-----------------------------------------+---------------------+---------------------+
|                                         |                     |                     |
+-----------------------------------------+---------------------+---------------------+

.. note::

   For this fix to work, the control variable must actually induce
   a change in a running LAMMPS simulation.  Typically this will only
   occur if there is some other command (e.g. a thermostat fix) which
   uses the control variable as an input parameter.  This could be done
   directly or indirectly, e.g. the other command uses a variable as
   input whose formula uses the control variable.  The other command
   should alter its behavior dynamically as the variable changes.

.. note::

   If there is a command you think could be used in this fashion,
   but does not currently allow a variable as an input parameter, please
   notify the LAMMPS developers.  It is often not difficult to enable a
   command to use a variable as an input parameter.

The group specified with this command is ignored.  However, note that
the process variable may be defined by calculations performed by
computes and fixes which store their own "group" definitions.

The PID controller is invoked once each *Nevery* timesteps.

The PID controller is implemented as a discretized version of
the following dynamic equation:

.. math::

   \frac{dc}{dt}  = \hat{E} -\alpha (K_p e + K_i \int_0^t e \, dt + K_d \frac{de}{dt} )

where *c* is the continuous time analog of the control variable,
*e* =\ *pvar*\ -\ *setpoint* is the error in the process variable, and
:math:`\alpha`, :math:`K_p`, :math:`K_i` , and :math:`K_d` are constants
set by the corresponding
keywords described above. The discretized version of this equation is:

.. math::

   c_n  = \hat{E} c_{n-1} -\alpha \left( K_p \tau e_n + K_i \tau^2 \sum_{i=1}^n e_i + K_d (e_n - e_{n-1}) \right)

where :math:`\tau = \mathtt{Nevery} \cdot \mathtt{timestep}` is the time
interval between updates,
and the subscripted variables indicate the values of *c* and *e* at
successive updates.

From the first equation, it is clear that if the three gain values
:math:`K_p`, :math:`K_i`, :math:`K_d` are dimensionless constants,
then :math:`\alpha` must have
units of [unit *cvar*\ ]/[unit *pvar*\ ]/[unit time] e.g. [ eV/K/ps
]. The advantage of this unit scheme is that the value of the
constants should be invariant under a change of either the MD timestep
size or the value of *Nevery*\ . Similarly, if the LAMMPS :doc:`unit style <units>` is changed, it should only be necessary to change
the value of :math:`\alpha` to reflect this, while leaving :math:`K_p`,
:math:`K_i`, and :math:`K_d` unaltered.

When choosing the values of the four constants, it is best to first
pick a value and sign for :math:`\alpha` that is consistent with the
magnitudes and signs of *pvar* and *cvar*\ .  The magnitude of :math:`K_p`
should then be tested over a large positive range keeping :math:`K_i = K_d =0`.
A good value for :math:`K_p` will produce a fast response in *pvar*\ ,
without overshooting the *setpoint*\ .  For many applications, proportional
feedback is sufficient, and so :math:`K_i` = K_d =0` can be used. In cases
where there is a substantial lag time in the response of *pvar* to a change
in *cvar*\ , this can be counteracted by increasing :math:`K_d`. In situations
where *pvar* plateaus without reaching *setpoint*\ , this can be
counteracted by increasing :math:`K_i`.  In the language of Charles Dickens,
:math:`K_p` represents the error of the present, :math:`K_i` the error of
the past, and :math:`K_d` the error yet to come.

Because this fix updates *cvar*\ , but does not initialize its value,
the initial value is that assigned by the user in the input script via
the :doc:`internal-style variable <variable>` command.  This value is
used (by the other LAMMPS command that used the variable) until this
fix performs its first update of *cvar* after *Nevery* timesteps.  On
the first update, the value of the derivative term is set to zero,
because the value of :math:`e_n-1` is not yet defined.

----------

The process variable *pvar* can be specified as the output of a
:doc:`compute <compute>` or :doc:`fix <fix>` or the evaluation of a
:doc:`variable <variable>`.  In each case, the compute, fix, or variable
must produce a global quantity, not a per-atom or local quantity.

If *pvar* begins with "c\_", a compute ID must follow which has been
previously defined in the input script and which generates a global
scalar or vector.  See the individual :doc:`compute <compute>` doc page
for details.  If no bracketed integer is appended, the scalar
calculated by the compute is used.  If a bracketed integer is
appended, the Ith value of the vector calculated by the compute is
used.  Users can also write code for their own compute styles and :doc:`add them to LAMMPS <Modify>`.

If *pvar* begins with "f\_", a fix ID must follow which has been
previously defined in the input script and which generates a global
scalar or vector.  See the individual :doc:`fix <fix>` doc page for
details.  Note that some fixes only produce their values on certain
timesteps, which must be compatible with when fix controller
references the values, or else an error results.  If no bracketed integer
is appended, the scalar calculated by the fix is used.  If a bracketed
integer is appended, the Ith value of the vector calculated by the fix
is used.  Users can also write code for their own fix style and :doc:`add them to LAMMPS <Modify>`.

If *pvar* begins with "v\_", a variable name must follow which has been
previously defined in the input script.  Only equal-style variables
can be referenced.  See the :doc:`variable <variable>` command for
details.  Note that variables of style *equal* define a formula which
can reference individual atom properties or thermodynamic keywords, or
they can invoke other computes, fixes, or variables when they are
evaluated, so this is a very general means of specifying the process
variable.

The target value *setpoint* for the process variable must be a numeric
value, in whatever units *pvar* is defined for.

The control variable *cvar* must be the name of an :doc:`internal-style variable <variable>` previously defined in the input script.  Note
that it is not specified with a "v\_" prefix, just the name of the
variable.  It must be an internal-style variable, because this fix
updates its value directly.  Note that other commands can use an
equal-style versus internal-style variable interchangeably.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix produces a global vector with 3 values which can be accessed
by various :doc:`output commands <Howto_output>`.  The values can be
accessed on any timestep, though they are only updated on timesteps
that are a multiple of *Nevery*\ .

The three values are the most recent updates made to the control
variable by each of the 3 terms in the PID equation above.  The first
value is the proportional term, the second is the integral term, the
third is the derivative term.

The units of the vector values will be whatever units the control
variable is in.  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix adapt <fix_adapt>`

**Default:** none
