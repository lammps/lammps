.. index:: fix ave/correlate

fix ave/correlate command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ave/correlate Nevery Nrepeat Nfreq value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/correlate = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of correlation time windows to accumulate
* Nfreq = calculate time window averages every this many timesteps
* one or more input values can be listed
* value = c_ID, c_ID[N], f_ID, f_ID[N], v_name

  .. parsed-literal::

       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID, I can include wildcard (see below)
       v_name = global value calculated by an equal-style variable with name
       v_name[I] = Ith component of a vector-style variable with name, I can include wildcard (see below)

* zero or more keyword/arg pairs may be appended
* keyword = *type* or *ave* or *start* or *prefactor* or *file* or *overwrite* or *title1* or *title2* or *title3*

  .. parsed-literal::

       *type* arg = *auto* or *upper* or *lower* or *auto/upper* or *auto/lower* or *full*
         auto = correlate each value with itself
         upper = correlate each value with each succeeding value
         lower = correlate each value with each preceding value
         auto/upper = auto + upper
         auto/lower = auto + lower
         full = correlate each value with every other value, including itself = auto + upper + lower
       *ave* args = *one* or *running*
         one = zero the correlation accumulation every Nfreq steps
         running = accumulate correlations continuously
       *start* args = Nstart
         Nstart = start accumulating correlations on this timestep
       *prefactor* args = value
         value = prefactor to scale all the correlation data by
       *file* arg = filename
         filename = name of file to output correlation data to
       *overwrite* arg = none = overwrite output file with only latest output
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *title3* arg = string
         string = text to print as 3rd line of output file

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/correlate 5 100 1000 c_myTemp file temp.correlate
   fix 1 all ave/correlate 1 50 10000 &
             c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] &
             type upper ave running title1 "My correlation data"

   fix 1 all ave/correlate 1 50 10000 c_thermo_press[*]

Description
"""""""""""

Use one or more global scalar values as inputs every few timesteps,
calculate time correlations between them at varying time intervals,
and average the correlation data over longer timescales.  The resulting
correlation values can be time integrated by
:doc:`variables <variable>` or used by other
:doc:`output commands <Howto_output>` such as
:doc:`thermo_style custom <thermo_style>`, and can also be written to a file.
See the :doc:`fix ave/correlate/long <fix_ave_correlate_long>` command for an
alternate method for computing correlation functions efficiently over
very long time windows.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an equal-style or vector-style
:doc:`variable <variable>`.  In each case, the compute, fix, or variable
must produce a global quantity, not a per-atom or local quantity.  If
you wish to spatial- or time-average or histogram per-atom quantities
from a compute, fix, or variable, then see the
:doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/atom <fix_ave_atom>`, or
:doc:`fix ave/histo <fix_ave_histo>` commands.  If you wish to convert a
per-atom quantity into a single global value, see the
:doc:`compute reduce <compute_reduce>` command.

The input values must be all scalars.  What kinds of
correlations between input values are calculated is determined by the
*type* keyword as discussed below.

:doc:`Computes <compute>` that produce global quantities are those which
do not have the word *atom* in their style name.  Only a few
:doc:`fixes <fix>` produce global quantities.  See the doc pages for
individual fixes for info on which ones produce such values.
:doc:`Variables <variable>` of style *equal* and *vector* are the only
ones that can be used with this fix.  Variables of style *atom* cannot
be used, since they produce per-atom values.

----------

For input values from a compute or fix or variable , the bracketed
index I can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the size of the vector, then an
asterisk with no numeric values means all indices from 1 to :math:`N`.  A
leading asterisk means all indices from 1 to n (inclusive).  A
trailing asterisk means all indices from m to :math:`N` (inclusive).
A middle asterisk means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual elements of the
vector had been listed one by one.  For example, the following two fix
ave/correlate commands are equivalent, since the :doc:`compute pressure
<compute_pressure>` command creates a global vector with six values:

.. code-block:: LAMMPS

   compute myPress all pressure NULL
   fix 1 all ave/correlate 1 50 10000 c_myPress[*]
   fix 1 all ave/correlate 1 50 10000 &
             c_myPress[1] c_myPress[2] c_myPress[3] &
             c_myPress[4] c_myPress[5] c_myPress[6]

.. note::

   For a vector-style variable, only the wildcard forms "\*n" or
   "m\*n" are allowed.  You must specify the upper bound, because
   vector-style variable lengths are not determined until the variable
   is evaluated.  If n is specified larger than the vector length
   turns out to be, zeroes are output for missing vector values.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and :math:`N_\text{freq}`
arguments specify on what timesteps the input values will be used to calculate
correlation data.  The input values are sampled every :math:`N_\text{every}`
time steps.  The correlation data for the preceding samples is computed on
time steps that are a multiple of :math:`N_\text{freq}`\ .  Consider a set of
samples from some initial time up to an output timestep.  The initial time
could be the beginning of the simulation or the last output time; see the *ave*
keyword for options.  For the set of samples, the correlation value
:math:`C_{ij}` is calculated as:

.. math::

   C_{ij}(\Delta t) = \left\langle V_i(t) V_j(t+\Delta t)\right\rangle,

which is the correlation value between input values :math:`V_i` and
:math:`V_j`, separated by time :math:`\Delta t`.  Note that the second value
:math:`V_j` in the pair is always the one sampled at the later time.  The
average is an average over every pair of samples in the set that are separated
by time :math:`\Delta t`.  The maximum :math:`\Delta t` used is of size
:math:`(N_\text{repeat} - 1) N_\text{every}`\ .
Thus the correlation between a pair of input values yields
:math:`N_\text{repeat}` correlation data:

.. math::

   C_{ij}(0), C_{ij}(N_\text{every}), C_{ij}(2N_\text{every}), \dotsc,
     C_{ij}\bigl((N_\text{repeat}-1) N_\text{every}\bigr)

For example, if :math:`N_\text{every}=5`, :math:`N_\text{repeat}=6`, and
:math:`N_\text{freq}=100`, then values on time steps
:math:`0, 5, 10, 15,\dotsc,100` will be used to compute the final averages
on time step 100.  Six averages will be computed: :math:`C_{ij}(0)`,
:math:`C_{ij}(5)`, :math:`C_{ij}(10)`, :math:`C_{ij}(15)`, :math:`C_{ij}(20)`,
and :math:`C_{ij}(25)`.  :math:`C_{ij}(10)` on time step 100 will
be the average of 19 samples, namely :math:`V_i(0) V_j(10)`,
:math:`V_i(5) V_j(15)`, :math:`V_i(10) V_j(20)`,
:math:`V_i(15) V_j(25), \dotsc,`
:math:`V_i(85) V_j(95)`, and :math:`V_i(90) V_j(100)`.

:math:`N_\text{freq}` must be a multiple of :math:`N_\text{every}`;
:math:`N_\text{every}` and :math:`N_\text{repeat}` must be non-zero.
Also, if the *ave* keyword is set to *one* which is the default, then
:math:`N_\text{freq} \ge (N_\text{repeat} -1) N_\text{every}` is required.

----------

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If no bracketed term is
appended, the global scalar calculated by the compute is used.  If a
bracketed term is appended, the :math:`I^\text{th}` element of the global
vector calculated by the compute is used.  See the discussion above for how
:math:`I` can be specified with a wildcard asterisk to effectively specify
multiple values.

Note that there is a :doc:`compute reduce <compute_reduce>` command
that can sum per-atom quantities into a global scalar or vector which
can then be accessed by fix ave/correlate.  It can also be a compute defined
not in your input script, but by :doc:`thermodynamic output <thermo_style>`
or other fixes such as :doc:`fix nvt <fix_nh>`
or :doc:`fix temp/rescale <fix_temp_rescale>`.  See the doc pages for
these commands which give the IDs of these computes.  Users can also
write code for their own compute styles and :doc:`add them to LAMMPS <Modify>`.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If no bracketed term is
appended, the global scalar calculated by the fix is used.  If a
bracketed term is appended, the :math:`I^\text{th}` element of the global
vector calculated by the fix is used.  See the discussion above for how
:math:`I` can be specified with a wildcard asterisk to effectively specify
multiple values.

Note that some fixes only produce their values on certain timesteps,
which must be compatible with :math:`N_\text{every}`, else an error will
result.  Users can also write code for their own fix styles and
:doc:`add them to LAMMPS <Modify>`.

If a value begins with "v\_", a variable name must follow which has been
previously defined in the input script.  Only equal-style or vector-style
variables can be referenced; the latter requires a bracketed term to specify
the :math:`I^\text{th}` element of the vector calculated by the variable.
See the :doc:`variable <variable>` command for details. Note that variables of
style *equal* or *vector* define a formula which can reference individual atom
properties or thermodynamic keywords, or they can invoke other computes, fixes,
or variables when they are evaluated, so this is a very general means of
specifying quantities to time correlate.

----------

Additional optional keywords also affect the operation of this fix.

The *type* keyword determines which pairs of input values are
correlated with each other.  For :math:`N` input values :math:`V_i`,
with :math:`i \in \{1,\dotsc,N\}`, let the number of pairs be
:math:`N_\text{pair}`.  Note that the second value in the
pair, :math:`V_i(t) V_j(t+\Delta t)`, is always the one sampled at the later
time.

* If *type* is set to *auto* then each input value is correlated with
  itself (i.e., :math:`C_{ii} = V_i^2` for :math:`i \in \{1,\dotsc,N\}`,
  so :math:`N_\text{pair} = N`).
* If *type* is set to *upper* then each input value is correlated with every
  succeeding value (i.e., :math:`C_{ij} = V_i V_j` for :math:`i < j`, so
  :math:`N_\text{pair} = N (N-1)/2`).
* If *type* is set to *lower* then each input value is correlated with every
  preceding value (i.e., :math:`C_{ij} = V_i V_j` for :math:`i > j`, so
  :math:`N_\text{pair} = N(N-1)/2`).
* If *type* is set to *auto/upper* then each input value is correlated
  with itself and every succeeding value (i.e., :math:`C_{ij} = V_i V_j`
  for :math:`i \ge j`, so :math:`N_\text{pair} = N(N+1)/2`).
* If *type* is set to *auto/lower* then each input value is correlated
  with itself and every preceding value (i.e., :math:`C_{ij} = V_i V_j`
  for :math:`i \le j`, so :math:`N_\text{pair} = N(N+1)/2`).
* If *type* is set to *full* then each input value is correlated with
  itself and every other value (i.e., :math:`C_{ij} = V_i V_j` for
  :math:`\{i,j\} = \{1,N\}`, so :math:`N_\text{pair} = N^2`).

The *ave* keyword determines what happens to the accumulation of correlation
samples every :math:`N_\text{freq}` timesteps.  If the *ave* setting is *one*,
then the accumulation is restarted or zeroed every :math:`N_\text{freq}`
timesteps.  Thus the outputs on successive :math:`N_\text{freq}` timesteps are
essentially independent of each other.  The exception is that the
:math:`C_{ij}(0) = V_i(t) V_j(t)` value at a time step :math:`t,` where
:math:`t` is a multiple of :math:`N_\text{freq}`, contributes to the
correlation output both at time :math:`t` and at time :math:`t+N_\text{freq}`.

If the *ave* setting is *running*, then the accumulation is never zeroed.
Thus the output of correlation data at any timestep is the average over samples
accumulated every :math:`N_\text{every}` steps since the fix was defined.
It can only be restarted by deleting the fix via the :doc:`unfix <unfix>`
command, or by re-defining the fix by re-specifying it.

The *start* keyword specifies what time step the accumulation of
correlation samples will begin on.  The default is step 0.  Setting it
to a larger value can avoid adding non-equilibrated data to the
correlation averages.

The *prefactor* keyword specifies a constant which will be used as a multiplier
on the correlation data after it is averaged.  It is effectively a scale factor
on :math:`V_i V_j`, which can be used to account for the size of the time
window or other unit conversions.

The *file* keyword allows a filename to be specified.  Every
:math:`N_\text{freq}` steps, an array of correlation data is written to the
file.  The number of rows is :math:`N_\text{repeat}`, as described above.
The number of columns is :math:`N_\text{pair}+2`, also as described above.
Thus the file ends up to be a series of these array sections.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one timestep worth of
output.  This option can only be used with the *ave running* setting.

The *title1*, *title2*, and *title3* keywords allow specification of
the strings that will be printed as the first three lines of the output file,
assuming the *file* keyword was used.  LAMMPS uses default values for each of
these, so they do not need to be specified.

By default, these header lines are as follows:

.. parsed-literal::

   # Time-correlated data for fix ID
   # TimeStep Number-of-time-windows
   # Index TimeDelta Ncount valueI\*valueJ valueI\*valueJ ...

In the first line, ID is replaced with the fix-ID.  The second line
describes the two values that are printed at the first of each section
of output.  In the third line the value pairs are replaced with the
appropriate fields from the fix ave/correlate command.

----------

Let :math:`S_{ij}` be a set of time correlation data for input values
:math:`I` and :math:`J`, namely the :math:`N_\text{repeat}` values:

.. math::

   S_{ij} = C_{ij}(0), C_{ij}(N_\text{every}), C_{ij}(2N_\text{every}),
    \dotsc, C_{ijI}\bigl((N_\text{repeat}-1) N_\text{every}\bigr)

As explained below, these data are output as one column of a global
array, which is effectively the correlation matrix.

The *trap* function defined for :doc:`equal-style variables <variable>`
can be used to perform a time integration of this vector of data,
using a trapezoidal rule.  This is useful for calculating various
quantities which can be derived from time correlation data.  If a
normalization factor is needed for the time integration, it can be
included in the variable formula or via the *prefactor* keyword.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to
:doc:`binary restart files <restart>`.  None of the
:doc:`fix_modify <fix_modify>` options are relevant to this fix.

This fix computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.  The values can only be
accessed on timesteps that are multiples of :math:`N_\text{freq}` since that is
when averaging is performed.  The global array has # of rows
:math:`N_\text{repeat}` and # of columns :math:`N_\text{pair}+2`.  The first
column has the time :math:`\Delta t` (in time steps) between the pairs of input
values used to calculate the correlation, as described above.  The second
column has the number of samples contributing to the correlation average, as
described above.  The remaining Npair columns are for :math:`I,J` pairs of the
:math:`N` input values, as determined by the *type* keyword, as described
above.

* For *type* = *auto*, the :math:`N_\text{pair} = N` columns are ordered:
  :math:`C_{11}, C_{22}, \dotsc, C_{NN}`
* For *type* = *upper*, the :math:`N_\text{pair} = N(N-1)/2` columns are
  ordered: :math:`C_{12}, C_{13}, \dotsc, C_{1N}, C_{23}, \dotsc, C_{2N},
  C_{34}, \dotsc, C_{N-1,N}`
* For *type* = *lower*, the :math:`N_\text{pair} = N(N-1)/2` columns are
  ordered: :math:`C_{21}, C_{31}, C_{32}, C_{41}, C_{42}, C_{43I}, \dotsc,
  C_{N1}, C_{N2}, \dotsc, C_{N,N-1}`
* For *type* = *auto/upper*, the :math:`N_\text{pair} = N(N+1)/2` columns are
  ordered: :math:`C_{11}, C_{12}, C_{13}, \dotsc, C_{1N}, C_{22}, C_{23},
  \dotsc, C_{2N}, C_{33}, C_{34}, \dotsc, C_{N-1,N}, C_{NN}`
* For *type* = *auto/lower*, the :math:`N_\text{pair} = N(N+1)/2` columns are
  ordered: :math:`C_{11}, C_{21}, C_{22}, C_{31}, C_{32}, C_{33}, C_{41},
  \dotsc, C_{44}, C_{N1}, C_{N2}, \dotsc, C_{N,N-1}, C_{NN}`
* For *type* = *full*, the :math:`N_\text{pair} = N^2` columns are ordered:
  :math:`C_{11}, C_{12}, \dotsc, C_{1N}, C_{21}, C_{22}, \dotsc, C_{2N},
  C_{31}, \dotsc, C_{3N}, \dotsc, C_{N1}, \dotsc, C_{N,N-1}, C_{NN}`

The array values calculated by this fix are treated as extensive.  If
you need to divide them by the number of atoms, you must do this in a
later processing step (e.g., when using them in a :doc:`variable <variable>`).

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix ave/correlate/long <fix_ave_correlate_long>`,
:doc:`compute <compute>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`fix ave/atom <fix_ave_atom>`, :doc:`fix ave/chunk <fix_ave_chunk>`,
:doc:`fix ave/histo <fix_ave_histo>`, :doc:`variable <variable>`

Default
"""""""

none

The option defaults are ave = one, type = auto, start = 0, no file
output, title 1,2,3 = strings as described above, and prefactor = 1.0.
