.. index:: fix pair/tracker

fix pair/tracker command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID pair/tracker N attribute1 attribute2 ... keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* pair/tracker = style name of this fix command
* N = prepare data for output every this many timesteps
* one or more attributes may be appended

  .. parsed-literal::

       possible attributes = id1 id2 time/created time/broken time/total
                             rmin rave x y z

  .. parsed-literal::

          id1, id2 = IDs of 2 atoms in each pair interaction
          time/created = the time the 2 atoms began interacting
          time/broken = the time the 2 atoms stopped interacting
          time/total = the total time the 2 atoms interacted
          r/min = the minimum radial distance between the 2 atoms during the interaction
          r/ave = the average radial distance between the 2 atoms during the interaction
          x, y, z = the center of mass position of the 2 atoms when they stopped interacting

* zero or more keyword/value pairs may be appended
* keyword = *time/min* or *type/include*

  .. parsed-literal::

       *time/min* value = T
         T = minimum interaction time 
       *type/include* value = arg1 arg2
         arg = separate lists of types (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all pair/tracker 1000 id1 id2 tmin 100
   fix 1 all pair/tracker 1000 time/created time/broken type 1 * type 2 3,4

Description
"""""""""""

Tracks properties of pairwise interactions between two atoms and records data
whenever the atoms move beyond the interaction cutoff or the interaction breaks. 
Must be used in conjuction with :doc:`pair tracker <pair_tracker>`.
Data is accumulated over a span of *N* timesteps after which it is cleared
The number of datums generated, aggregated across all processors, equals 
the number of broken interactions. Interactions are only included if both
atoms are in the included in the specified fix group. Additional filters can be
applied using the *tmin* or *type* keywords described below.

.. note::

   For extremely long-lived interactions, the calculation of *r/ave* may not be 
   correct due to double overflow.

The *time/min* keyword defines a minimum amount of time atoms have to be interacting
to save data. This can be used to censor short-lived interactions. The *type/include*
keyword filters interactions based on the types of the two atoms. Data is 
only saved for interactions between atoms with types in the two lists. 
Each list consists of a series of type
ranges separated by commas. The range can be specified as a
single numeric value, or a wildcard asterisk can be used to specify a range
of values.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  For
example, if M = the number of atom types, then an asterisk with no numeric
values means all types from 1 to M.  A leading asterisk means all types
from 1 to n (inclusive).  A trailing asterisk means all types from n to M
(inclusive).  A middle asterisk means all types from m to n (inclusive).
Note that all atom types must be included in exactly one of the N collections.
Multiple *type/include* keywords may be added.

----------

Restart, fix_modify, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  
None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Output info
"""""""""""

This compute calculates a local vector or local array depending on the
number of input values.  The length of the vector or number of rows in
the array is the number of bonds, angles, etc.  If a single input is
specified, a local vector is produced.  If two or more inputs are
specified, a local array is produced where the number of columns = the
number of inputs.  The vector or array can be accessed by any command
that uses local values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The vector or array values will be doubles that correspond to the
specified attribute.

Restrictions
""""""""""""
 none
 
Related commands
""""""""""""""""

:doc:`pair tracker <pair_tracker>`

Default
"""""""

none
