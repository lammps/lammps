.. index:: compute global/atom

compute global/atom command
===========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID style index input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* global/atom = style name of this compute command
* index = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       c_ID = per-atom vector calculated by a compute with ID
       c_ID[I] = Ith column of per-atom array calculated by a compute with ID
       f_ID = per-atom vector calculated by a fix with ID
       f_ID[I] = Ith column of per-atom array calculated by a fix with ID
       v_name = per-atom vector calculated by an atom-style variable with name

* one or more inputs can be listed
* input = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       c_ID = global vector calculated by a compute with ID
       c_ID[I] = Ith column of global array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global vector calculated by a fix with ID
       f_ID[I] = Ith column of global array calculated by a fix with ID, I can include wildcard (see below)
       v_name = global vector calculated by a vector-style variable with name



Examples
""""""""


.. parsed-literal::

   compute 1 all global/atom c_chunk c_com[1\] c_com[2\] c_com[3\]
   compute 1 all global/atom c_chunk c_com[\*\]

Description
"""""""""""

Define a calculation that assigns global values to each atom from
vectors or arrays of global values.  The specified *index* parameter
is used to determine which global value is assigned to each atom.

The *index* parameter must reference a per-atom vector or array from a
:doc:`compute <compute>` or :doc:`fix <fix>` or the evaluation of an
atom-style :doc:`variable <variable>`.  Each *input* value must
reference a global vector or array from a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an vector-style
:doc:`variable <variable>`.  Details are given below.

The *index* value for an atom is used as a index I (from 1 to N) into
the vector associated with each of the input values.  The Ith value
from the input vector becomes one output value for that atom.  If the
atom is not in the specified group, or the index I < 1 or I > M, where
M is the actual length of the input vector, then an output value of
0.0 is assigned to the atom.

An example of how this command is useful, is in the context of
"chunks" which are static or dynamic subsets of atoms.  The :doc:`compute chunk/atom <compute_chunk_atom>` command assigns unique chunk IDs
to each atom.  It's output can be used as the *index* parameter for
this command.  Various other computes with "chunk" in their style
name, such as :doc:`compute com/chunk <compute_com_chunk>` or :doc:`compute msd/chunk <compute_msd_chunk>`, calculate properties for each
chunk.  The output of these commands are global vectors or arrays,
with one or more values per chunk, and can be used as input values for
this command.  This command will then assign the global chunk value to
each atom in the chunk, producing a per-atom vector or per-atom array
as output.  The per-atom values can then be output to a dump file or
used by any command that uses per-atom values from a compute as input,
as discussed on the :doc:`Howto output <Howto_output>` doc page.

As a concrete example, these commands will calculate the displacement
of each atom from the center-of-mass of the molecule it is in, and
dump those values to a dump file.  In this case, each molecule is a
chunk.


.. parsed-literal::

   compute cc1 all chunk/atom molecule
   compute myChunk all com/chunk cc1
   compute prop all property/atom xu yu zu
   compute glob all global/atom c_cc1 c_myChunk[\*]
   variable dx atom c_prop[1]-c_glob[1]
   variable dy atom c_prop[2]-c_glob[2]
   variable dz atom c_prop[3]-c_glob[3]
   variable dist atom sqrt(v_dx\*v_dx+v_dy\*v_dy+v_dz\*v_dz)
   dump 1 all custom 100 tmp.dump id xu yu zu c_glob[1] c_glob[2] c_glob[3] &
        v_dx v_dy v_dz v_dist
   dump_modify 1 sort id

You can add these commands to the bench/in.chain script to see how
they work.


----------


Note that for input values from a compute or fix, the bracketed index
I can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or "\*n"
or "n\*" or "m\*n".  If N = the size of the vector (for *mode* = scalar)
or the number of columns in the array (for *mode* = vector), then an
asterisk with no numeric values means all indices from 1 to N.  A
leading asterisk means all indices from 1 to n (inclusive).  A
trailing asterisk means all indices from n to N (inclusive).  A middle
asterisk means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 compute global/atom commands
are equivalent, since the :doc:`compute com/chunk <compute_com_chunk>`
command creates a global array with 3 columns:


.. parsed-literal::

   compute cc1 all chunk/atom molecule
   compute com all com/chunk cc1
   compute 1 all global/atom c_cc1 c_com[1] c_com[2] c_com[3]
   compute 1 all global/atom c_cc1 c_com[\*]


----------


This section explains the *index* parameter.  Note that it must
reference per-atom values, as contrasted with the *input* values which
must reference global values.

Note that all of these options generate floating point values.  When
they are used as an index into the specified input vectors, they
simple rounded down to convert the value to integer indices.  The
final values should range from 1 to N (inclusive), since they are used
to access values from N-length vectors.

If *index* begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  The compute must generate
per-atom quantities.  See the individual :doc:`compute <compute>` doc
page for details.  If no bracketed integer is appended, the per-atom
vector calculated by the compute is used.  If a bracketed integer is
appended, the Ith column of the per-atom array calculated by the
compute is used.  Users can also write code for their own compute
styles and :doc:`add them to LAMMPS <Modify>`.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

If *index* begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  The Fix must generate
per-atom quantities.  See the individual :doc:`fix <fix>` doc page for
details.  Note that some fixes only produce their values on certain
timesteps, which must be compatible with when compute global/atom
references the values, else an error results.  If no bracketed integer
is appended, the per-atom vector calculated by the fix is used.  If a
bracketed integer is appended, the Ith column of the per-atom array
calculated by the fix is used.  Users can also write code for their
own fix style and :doc:`add them to LAMMPS <Modify>`.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

If *index* begins with "v\_", a variable name must follow which has
been previously defined in the input script.  It must be an
:doc:`atom-style variable <variable>`.  Atom-style variables can
reference thermodynamic keywords and various per-atom attributes, or
invoke other computes, fixes, or variables when they are evaluated, so
this is a very general means of generating per-atom quantities to use
as *index*\ .


----------


This section explains the kinds of *input* values that can be used.
Note that inputs reference global values, as contrasted with the
*index* parameter which must reference per-atom values.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  The compute must generate a
global vector or array.  See the individual :doc:`compute <compute>` doc
page for details.  If no bracketed integer is appended, the vector
calculated by the compute is used.  If a bracketed integer is
appended, the Ith column of the array calculated by the compute is
used.  Users can also write code for their own compute styles and :doc:`add them to LAMMPS <Modify>`.  See the discussion above for how
I can be specified with a wildcard asterisk to effectively specify
multiple values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  The fix must generate a
global vector or array.  See the individual :doc:`fix <fix>` doc page
for details.  Note that some fixes only produce their values on
certain timesteps, which must be compatible with when compute
global/atom references the values, else an error results.  If no
bracketed integer is appended, the vector calculated by the fix is
used.  If a bracketed integer is appended, the Ith column of the array
calculated by the fix is used.  Users can also write code for their
own fix style and :doc:`add them to LAMMPS <Modify>`.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  It must be a
:doc:`vector-style variable <variable>`.  Vector-style variables can
reference thermodynamic keywords and various other attributes of
atoms, or invoke other computes, fixes, or variables when they are
evaluated, so this is a very general means of generating a vector of
global quantities which the *index* parameter will reference for
assignment of global values to atoms.


----------


**Output info:**

If a single input is specified this compute produces a per-atom
vector.  If multiple inputs are specified, this compute produces a
per-atom array values, where the number of columns is equal to the
number of inputs specified.  These values can be used by any command
that uses per-atom vector or array values from a compute as input.
See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector or array values will be in whatever units the
corresponding input values are in.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix <fix>`, :doc:`variable <variable>`,
:doc:`compute chunk/atom <compute_chunk_atom>`, :doc:`compute reduce <compute_reduce>`

**Default:** none
