.. index:: compute reduce

compute reduce command
======================

compute reduce/region command
=============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID style arg mode input1 input2 ... keyword args ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *reduce* or *reduce/region*
  
  .. parsed-literal::
  
       *reduce* arg = none
       *reduce/region* arg = region-ID
         region-ID = ID of region to use for choosing atoms

* mode = *sum* or *min* or *max* or *ave* or *sumsq* or *avesq*
* one or more inputs can be listed
* input = x, y, z, vx, vy, vz, fx, fy, fz, c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       x,y,z,vx,vy,vz,fx,fy,fz = atom attribute (position, velocity, force component)
       c_ID = per-atom or local vector calculated by a compute with ID
       c_ID[I] = Ith column of per-atom or local array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = per-atom or local vector calculated by a fix with ID
       f_ID[I] = Ith column of per-atom or local array calculated by a fix with ID, I can include wildcard (see below)
       v_name = per-atom vector calculated by an atom-style variable with name

* zero or more keyword/args pairs may be appended
* keyword = *replace*
  
  .. parsed-literal::
  
       *replace* args = vec1 vec2
         vec1 = reduced value from this input vector will be replaced
         vec2 = replace it with vec1[N] where N is index of max/min value from vec2



Examples
""""""""


.. parsed-literal::

   compute 1 all reduce sum c_force
   compute 1 all reduce/region subbox sum c_force
   compute 2 all reduce min c_press[2] f_ave v_myKE
   compute 2 all reduce min c_press[\*] f_ave v_myKE
   compute 3 fluid reduce max c_index[1] c_index[2] c_dist replace 1 3 replace 2 3

Description
"""""""""""

Define a calculation that "reduces" one or more vector inputs into
scalar values, one per listed input.  The inputs can be per-atom or
local quantities; they cannot be global quantities.  Atom attributes
are per-atom quantities, :doc:`computes <compute>` and :doc:`fixes <fix>`
may generate any of the three kinds of quantities, and :doc:`atom-style variables <variable>` generate per-atom quantities.  See the
:doc:`variable <variable>` command and its special functions which can
perform the same operations as the compute reduce command on global
vectors.

The reduction operation is specified by the *mode* setting.  The *sum*
option adds the values in the vector into a global total.  The *min*
or *max* options find the minimum or maximum value across all vector
values.  The *ave* setting adds the vector values into a global total,
then divides by the number of values in the vector.  The *sumsq*
option sums the square of the values in the vector into a global
total.  The *avesq* setting does the same as *sumsq*\ , then divides the
sum of squares by the number of values.  The last two options can be
useful for calculating the variance of some quantity, e.g. variance =
sumsq - ave\^2.

Each listed input is operated on independently.  For per-atom inputs,
the group specified with this command means only atoms within the
group contribute to the result.  For per-atom inputs, if the compute
reduce/region command is used, the atoms must also currently be within
the region.  Note that an input that produces per-atom quantities may
define its own group which affects the quantities it returns.  For
example, if a compute is used as an input which generates a per-atom
vector, it will generate values of 0.0 for atoms that are not in the
group specified for that compute.

Each listed input can be an atom attribute (position, velocity, force
component) or can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an atom-style
:doc:`variable <variable>`.

Note that for values from a compute or fix, the bracketed index I can
be specified using a wildcard asterisk with the index to effectively
specify multiple values.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the size of the vector (for *mode* = scalar) or the
number of columns in the array (for *mode* = vector), then an asterisk
with no numeric values means all indices from 1 to N.  A leading
asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to N (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 compute reduce commands are
equivalent, since the :doc:`compute stress/atom <compute_stress_atom>`
command creates a per-atom array with 6 columns:


.. parsed-literal::

   compute myPress all stress/atom NULL
   compute 2 all reduce min c_myPress[\*]
   compute 2 all reduce min c_myPress[1] c_myPress[2] c_myPress[3] &
                            c_myPress[4] c_myPress[5] c_myPress[6]


----------


The atom attribute values (x,y,z,vx,vy,vz,fx,fy,fz) are
self-explanatory.  Note that other atom attributes can be used as
inputs to this fix by using the :doc:`compute property/atom <compute_property_atom>` command and then specifying
an input value from that compute.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  Computes can generate
per-atom or local quantities.  See the individual
:doc:`compute <compute>` doc page for details.  If no bracketed integer
is appended, the vector calculated by the compute is used.  If a
bracketed integer is appended, the Ith column of the array calculated
by the compute is used.  Users can also write code for their own
compute styles and :doc:`add them to LAMMPS <Modify>`.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  Fixes can generate per-atom
or local quantities.  See the individual :doc:`fix <fix>` doc page for
details.  Note that some fixes only produce their values on certain
timesteps, which must be compatible with when compute reduce
references the values, else an error results.  If no bracketed integer
is appended, the vector calculated by the fix is used.  If a bracketed
integer is appended, the Ith column of the array calculated by the fix
is used.  Users can also write code for their own fix style and :doc:`add them to LAMMPS <Modify>`.  See the discussion above for how I can
be specified with a wildcard asterisk to effectively specify multiple
values.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  It must be an
:doc:`atom-style variable <variable>`.  Atom-style variables can
reference thermodynamic keywords and various per-atom attributes, or
invoke other computes, fixes, or variables when they are evaluated, so
this is a very general means of generating per-atom quantities to
reduce.


----------


If the *replace* keyword is used, two indices *vec1* and *vec2* are
specified, where each index ranges from 1 to the # of input values.
The replace keyword can only be used if the *mode* is *min* or *max*\ .
It works as follows.  A min/max is computed as usual on the *vec2*
input vector.  The index N of that value within *vec2* is also stored.
Then, instead of performing a min/max on the *vec1* input vector, the
stored index is used to select the Nth element of the *vec1* vector.

Thus, for example, if you wish to use this compute to find the bond
with maximum stretch, you can do it as follows:


.. parsed-literal::

   compute 1 all property/local batom1 batom2
   compute 2 all bond/local dist
   compute 3 all reduce max c_1[1] c_1[2] c_2 replace 1 3 replace 2 3
   thermo_style custom step temp c_3[1] c_3[2] c_3[3]

The first two input values in the compute reduce command are vectors
with the IDs of the 2 atoms in each bond, using the :doc:`compute property/local <compute_property_local>` command.  The last input
value is bond distance, using the :doc:`compute bond/local <compute_bond_local>` command.  Instead of taking the
max of the two atom ID vectors, which does not yield useful
information in this context, the *replace* keywords will extract the
atom IDs for the two atoms in the bond of maximum stretch.  These atom
IDs and the bond stretch will be printed with thermodynamic output.


----------


If a single input is specified this compute produces a global scalar
value.  If multiple inputs are specified, this compute produces a
global vector of values, the length of which is equal to the number of
inputs specified.

As discussed below, for the *sum* and *sumsq* modes, the value(s)
produced by this compute are all "extensive", meaning their value
scales linearly with the number of atoms involved.  If normalized
values are desired, this compute can be accessed by the :doc:`thermo_style custom <thermo_style>` command with :doc:`thermo_modify norm yes <thermo_modify>` set as an option.  Or it can be accessed by a
:doc:`variable <variable>` that divides by the appropriate atom count.


----------


**Output info:**

This compute calculates a global scalar if a single input value is
specified or a global vector of length N where N is the number of
inputs, and which can be accessed by indices 1 to N.  These values can
be used by any command that uses global scalar or vector values from a
compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

All the scalar or vector values calculated by this compute are
"intensive", except when the *sum* or *sumsq* modes are used on
per-atom or local vectors, in which case the calculated values are
"extensive".

The scalar or vector values will be in whatever :doc:`units <units>` the
quantities being reduced are in.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix <fix>`, :doc:`variable <variable>`

**Default:** none


