.. index:: fix external

fix external command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID external mode args

* ID, group-ID are documented in :doc:`fix <fix>` command
* external = style name of this fix command
* mode = *pf/callback* or *pf/array*
  
  .. parsed-literal::
  
       *pf/callback* args = Ncall Napply
         Ncall = make callback every Ncall steps
         Napply = apply callback forces every Napply steps
       *pf/array* args = Napply
         Napply = apply array forces every Napply steps



Examples
""""""""


.. parsed-literal::

   fix 1 all external pf/callback 1 1
   fix 1 all external pf/callback 100 1
   fix 1 all external pf/array 10

Description
"""""""""""

This fix allows external programs that are running LAMMPS through its
:doc:`library interface <Howto_library>` to modify certain LAMMPS
properties on specific timesteps, similar to the way other fixes do.
The external driver can be a :doc:`C/C++ or Fortran program <Howto_library>` or a :doc:`Python script <Python_head>`.


----------


If mode is *pf/callback* then the fix will make a callback every
*Ncall* timesteps or minimization iterations to the external program.
The external program computes forces on atoms by setting values in an
array owned by the fix.  The fix then adds these forces to each atom
in the group, once every *Napply* steps, similar to the way the :doc:`fix addforce <fix_addforce>` command works.  Note that if *Ncall* >
*Napply*\ , the force values produced by one callback will persist, and
be used multiple times to update atom forces.

The callback function "foo" is invoked by the fix as:


.. parsed-literal::

   foo(void \*ptr, bigint timestep, int nlocal, int \*ids, double \*\*x, double \*\*fexternal);

The arguments are as follows:

* ptr = pointer provided by and simply passed back to external driver
* timestep = current LAMMPS timestep
* nlocal = # of atoms on this processor
* ids = list of atom IDs on this processor
* x = coordinates of atoms on this processor
* fexternal = forces to add to atoms on this processor

Note that timestep is a "bigint" which is defined in src/lmptype.h,
typically as a 64-bit integer.

Fexternal are the forces returned by the driver program.

The fix has a set\_callback() method which the external driver can call
to pass a pointer to its foo() function.  See the
couple/lammps\_quest/lmpqst.cpp file in the LAMMPS distribution for an
example of how this is done.  This sample application performs
classical MD using quantum forces computed by a density functional
code `Quest <quest_>`_.

.. _quest: http://dft.sandia.gov/Quest




----------


If mode is *pf/array* then the fix simply stores force values in an
array.  The fix adds these forces to each atom in the group, once
every *Napply* steps, similar to the way the :doc:`fix addforce <fix_addforce>` command works.

The name of the public force array provided by the FixExternal
class is


.. parsed-literal::

   double \*\*fexternal;

It is allocated by the FixExternal class as an (N,3) array where N is
the number of atoms owned by a processor.  The 3 corresponds to the
fx, fy, fz components of force.

It is up to the external program to set the values in this array to
the desired quantities, as often as desired.  For example, the driver
program might perform an MD run in stages of 1000 timesteps each.  In
between calls to the LAMMPS :doc:`run <run>` command, it could retrieve
atom coordinates from LAMMPS, compute forces, set values in fexternal,
etc.


----------


To use this fix during energy minimization, the energy corresponding
to the added forces must also be set so as to be consistent with the
added forces.  Otherwise the minimization will not converge correctly.

This can be done from the external driver by calling this public
method of the FixExternal class:


.. parsed-literal::

   void set_energy(double eng);

where eng is the potential energy.  Eng is an extensive quantity,
meaning it should be the sum over per-atom energies of all affected
atoms.  It should also be provided in :doc:`energy units <units>`
consistent with the simulation.  See the details below for how to
insure this energy setting is used appropriately in a minimization.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the potential "energy" set by the external driver to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.  This is a fictitious quantity but is
needed so that the :doc:`minimize <minimize>` command can include the
forces added by this fix in a consistent manner.  I.e. there is a
decrease in potential energy when atoms move in the direction of the
added force.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the interactions computed by the
external program to the system's virial as part of :doc:`thermodynamic output <thermo_style>`. The default is *virial yes*

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy discussed above.  The scalar stored by this fix is "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the fictitious potential energy associated with the
   added forces to be included in the total potential energy of the
   system (the quantity being minimized), you MUST enable the
   :doc:`fix_modify <fix_modify>` *energy* option for this fix.

Restrictions
""""""""""""
 none

**Related commands:** none

**Default:** none


