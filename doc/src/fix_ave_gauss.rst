.. index:: fix ave/gauss

fix ave/gauss command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ave/gauss Nfreq Nwindow value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/gauss = style name of this fix command
* Nfreq = update statistics every this many time steps
* Nwindow = compute based on this many most recent steps
* one or more input variables can be listed
* value = v_name

  .. parsed-literal::

       v_name = value(s) calculated by an equal-style variable with name
       v_name[I] = value calculated by a vector-style variable with name

* zero or more keyword/arg pairs may be appended
* keyword = *delay* 

  .. parsed-literal::

       *delay* args = Ndelay
         Ndelay = keep the results from Ndelay updates ago (keyword can be repeated)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/gauss 100 1000 v_volume
   fix 1 all ave/gauss 1000 200 v_volume delay 10

Description
"""""""""""

Using one or more variables as input every few time steps, calculate the
average and standard deviation.  Calulations are performed similar to
:doc:`fix ave/time <fix_ave_time>` in windowed scalar mode.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value must be an equal-style or indexed vector-style
:doc:`variable <variable>`.  The variable must produce a global scalar
quantity, not a per-atom or local quantity.  Only a few :doc:`fixes <fix>`
produce global quantities.  See the doc pages for individual fixes for info on
which ones produce such values.  Variables of style *atom* cannot
be used, since they produce per-atom values.

The result of this fix can be accessed as a vector, containing the interleaved
average and standard deviation of each input in order.  The average of input 1
has index 1, the standard deviation index 2, the average of input 2 has index 3
and so on.

----------

The :math:`N_\text{freq}` and :math:`N_\text{window}` arguments specify which values
will be used in order to contribute to the average.  The final averaged quantities
are generated on time steps that are a multiple of :math:`N_\text{freq}`\ .
The average is over the value of the quantities computed at the most recent
:math:`N_\text{window}` time steps.  If the variable accesses some compute values,
these must be valid at all timesteps evaluated.

The values need not have any special relation: it is valid to have a window
larger than :math:`N_\text{freq}` as well as the other way around.
For example, if :math:`N_\text{freq}=100` and :math:`N_\text{window}=5`,
then values from time steps 96, 97, 98, 99, and 100 will be used.
This mean some intervening time  steps do not contribute to the result.
If :math:`N_\text{freq}=5` and :math:`N_\text{window}=10`, then values will
first be calculated on step 5 from steps 1-5, on step 10 from 1-10, on
step 15 from 5-15 and so on, forming a rolling average.

----------

Variable names must be given beginning with "v\_" and the variable name 
following.  The variable must have been previously defined in the input script.
Only equal-style or vector-style variables can be used, which both
produce global values.  Vector-style variables require a bracketed term
to specify the Ith element of the vector calculated by the variable.

Note that variables of style *equal* and *vector* define a formula
which can reference individual atom properties or thermodynamic
keywords, or they can invoke other computes, fixes, or variables when
they are evaluated, so this is a very general means of specifying
quantities to time average.

----------

Additional optional keywords also affect the operation of this fix.

The *delay* allows keeping a record of previous results of averaging.
Only if no *delay* keyword is given, *delay 0* is implied.

For example, this will output values which are delayed by 10 invocations,
meaning 10000 time steps:

.. code-block:: LAMMPS

   fix 1 all ave/gauss 1000 200 v_volume delay 10

For each instance of *delay* given, an additional row on the fix output
array can be accessed, containing the N-th last evaluation result.  These
indices are based on the order of the keywords given, not their numerical values.
For example, using "delay 0 delay 10 delay 5", the most recent average of the first
input value would be accessed as "f_name[1][1]", "f_name[1][2]" is the 10th most
recent and "f_name[1][3]" is the 5th most recent.  Vector access is always the
same as the first array row, corresponding to the first *delay* keyword.  If
the most recent result is desired, *delay 0* must be specified explicitly.

This fix can be used in conjunction with :doc:`fix halt <fix_halt>` to stop
a run automatically if a quantity is converged to within some limit:

.. code-block:: LAMMPS

   variable target equal etot
   fix aveg all ave/gauss 1000 200 v_target delay 0 delay 10
   variable stopcond equal "abs(f_aveg[1]-f_aveg[1][2])<f_aveg[2]"
   fix fhalt all halt 1000 v_stopcond == 1

In this example, every 1000 time steps, the average and standard deviation
of the total energy over the previous 200 time steps are calculated.  If the
difference between the most recent and 10-th most recent average is lower than
the most recent standard deviation, the run is stopped.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

This fix produces a global vector and global array which can be accessed
by various :doc:`output commands <Howto_output>`.
The values can be accessed on any time step, but may not be current.

A vector is produced with # of elements = 2 * number of inputts.  Odd-numbered
elements contain the average, even-numbered elements contain the standard
deviation of the value.  An array is produced having # of rows = number
of *delay* keywords given and # of columns = 2 * number of inputs, using
the same ordering as the vector output.

Each element can be either "intensive" or "extensive", depending on whether
the values contributing to the element are "intensive" or "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This compute is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix ave/time <fix_ave_time>`,

Default
"""""""

The option defaults are delay 0.
