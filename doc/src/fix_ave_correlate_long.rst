.. index:: fix ave/correlate/long

fix ave/correlate/long command
==============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID ave/correlate/long Nevery Nfreq value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/correlate/long = style name of this fix command
* Nevery = use input values every this many timesteps
* Nfreq = save state of the time correlation functions every this many timesteps
* one or more input values can be listed
* value = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID
       v_name = global value calculated by an equal-style variable with name

* zero or more keyword/arg pairs may be appended
* keyword = *type* or *start* or *file* or *overwrite* or *title1* or *title2* or *ncorr* or *p* or *m*
  
  .. parsed-literal::
  
       *type* arg = *auto* or *upper* or *lower* or *auto/upper* or *auto/lower* or *full*
         auto = correlate each value with itself
         upper = correlate each value with each succeeding value
         lower = correlate each value with each preceding value
         auto/upper = auto + upper
         auto/lower = auto + lower
         full = correlate each value with every other value, including itself = auto + upper + lower
       *start* args = Nstart
         Nstart = start accumulating correlations on this timestep
       *file* arg = filename
         filename = name of file to output correlation data to
       *overwrite* arg = none = overwrite output file with only latest output
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *ncorr* arg = Ncorrelators
         Ncorrelators = number of correlators to store
       *nlen* args = Nlen
         Nlen = length of each correlator
       *ncount* args = Ncount
         Ncount = number of values over which succesive correlators are averaged



Examples
""""""""


.. parsed-literal::

   fix 1 all ave/correlate/long 5 1000 c_myTemp file temp.correlate
   fix 1 all ave/correlate/long 1 10000 &
             c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] &
             type upper title1 "My correlation data" nlen 15 ncount 3

Description
"""""""""""

This fix is similar in spirit and syntax to the :doc:`fix ave/correlate <fix_ave_correlate>`.  However, this fix allows the
efficient calculation of time correlation functions on the fly over
extremely long time windows without too much CPU overhead, using a
multiple-tau method :ref:`(Ramirez) <Ramirez>` that decreases the resolution
of the stored correlation function with time.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value can be the result of a compute or fix or the
evaluation of an equal-style variable. See the :doc:`fix ave/correlate <fix_ave_correlate>` doc page for details.

The *Nevery* and *Nfreq* arguments specify on what timesteps the input
values will be used to calculate correlation data, and the frequency
with which the time correlation functions will be output to a file.
Note that there is no *Nrepeat* argument, unlike the :doc:`fix ave/correlate <fix_ave_correlate>` command.

The optional keywords *ncorr*\ , *nlen*\ , and *ncount* are unique to this
command and determine the number of correlation points calculated and
the memory and CPU overhead used by this calculation. *Nlen* and
*ncount* determine the amount of averaging done at longer correlation
times.  The default values *nlen=16*\ , *ncount=2* ensure that the
systematic error of the multiple-tau correlator is always below the
level of the statistical error of a typical simulation (which depends
on the ensemble size and the simulation length).

The maximum correlation time (in time steps) that can be reached is
given by the formula (nlen-1) \* ncount\^(ncorr-1).  Longer correlation
times are discarded and not calculated.  With the default values of
the parameters (ncorr=20, nlen=16 and ncount=2), this corresponds to
7864320 time steps.  If longer correlation times are needed, the value
of ncorr should be increased. Using nlen=16 and ncount=2, with
ncorr=30, the maximum number of steps that can be correlated is
80530636808.  If ncorr=40, correlation times in excess of 8e12 time
steps can be calculated.

The total memory needed for each correlation pair is roughly
4\*ncorr\*nlen\*8 bytes. With the default values of the parameters, this
corresponds to about 10 KB.

For the meaning of the additional optional keywords, see the :doc:`fix ave/correlate <fix_ave_correlate>` doc page.

**Restart, fix\_modify, output, run start/stop, minimize info:**

Since this fix in intended for the calculation of time correlation
functions over very long MD simulations, the information about this
fix is written automatically to binary restart files, so that the time
correlation calculation can continue in subsequent simulations. None
of the fix\_modify options are relevant to this fix.

No parameter of this fix can be used with the start/stop keywords of
the run command. This fix is not invoked during energy minimization.

Restrictions
""""""""""""


This compute is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix ave/correlate <fix_ave_correlate>`

**Default:** none

The option defaults for keywords that are also keywords for the :doc:`fix ave/correlate <fix_ave_correlate>` command are as follows: type =
auto, start = 0, no file output, title 1,2 = strings as described on
the :doc:`fix ave/correlate <fix_ave_correlate>` doc page.

The option defaults for keywords unique to this command are as
follows: ncorr=20, nlen=16, ncount=2.


----------


.. _Ramirez:



**(Ramirez)** J. Ramirez, S.K. Sukumaran, B. Vorselaars and
A.E. Likhtman, J. Chem. Phys. 133, 154103 (2010).


