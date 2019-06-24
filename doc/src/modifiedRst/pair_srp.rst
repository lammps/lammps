.. index:: pair\_style srp

pair\_style srp command
=======================

Syntax
""""""

pair\_style srp cutoff btype dist keyword value ...

* cutoff = global cutoff for SRP interactions (distance units)
* btype = bond type to apply SRP interactions to (can be wildcard, see below)
* distance = *min* or *mid*
* zero or more keyword/value pairs may be appended
* keyword = *exclude*
  
  .. parsed-literal::
  
       *bptype* value = atom type for bond particles
       *exclude* value = *yes* or *no*



Examples
""""""""


.. parsed-literal::

   pair_style hybrid dpd 1.0 1.0 12345 srp 0.8 1 mid exclude yes
   pair_coeff 1 1 dpd 60.0 4.5 1.0
   pair_coeff 1 2 none
   pair_coeff 2 2 srp 100.0 0.8

   pair_style hybrid dpd 1.0 1.0 12345 srp 0.8 \* min exclude yes
   pair_coeff 1 1 dpd 60.0 50 1.0
   pair_coeff 1 2 none
   pair_coeff 2 2 srp 40.0

   pair_style hybrid srp 0.8 2 mid
   pair_coeff 1 1 none
   pair_coeff 1 2 none
   pair_coeff 2 2 srp 100.0 0.8

Description
"""""""""""

Style *srp* computes a soft segmental repulsive potential (SRP) that
acts between pairs of bonds. This potential is useful for preventing
bonds from passing through one another when a soft non-bonded
potential acts between beads in, for example, DPD polymer chains.  An
example input script that uses this command is provided in
examples/USER/srp.

Bonds of specified type *btype* interact with one another through a
bond-pairwise potential, such that the force on bond *i* due to bond
*j* is as follows

.. math::

   F^{SRP}_{ij} & = & C(1-r/r_c)\hat{r}_{ij} \qquad r < r_c \\


where *r* and *rij* are the distance and unit vector between the two
bonds.  Note that *btype* can be specified as an asterisk "\*", which
case the interaction is applied to all bond types.  The *mid* option
computes *r* and *rij* from the midpoint distance between bonds. The
*min* option computes *r* and *rij* from the minimum distance between
bonds. The force acting on a bond is mapped onto the two bond atoms
according to the lever rule,

.. math::

   F_{i1}^{SRP} & = & F^{SRP}_{ij}(L) \\
   F_{i2}^{SRP} & = & F^{SRP}_{ij}(1-L)  


where *L* is the normalized distance from the atom to the point of
closest approach of bond *i* and *j*\ . The *mid* option takes *L* as
0.5 for each interaction as described in :ref:`(Sirk) <Sirk2>`.

The following coefficients must be defined via the
:doc:`pair\_coeff <pair_coeff>` command as in the examples above, or in
the data file or restart file read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* *C* (force units)
* *rc* (distance units)

The last coefficient is optional. If not specified, the global cutoff
is used.

.. note::

   Pair style srp considers each bond of type *btype* to be a
   fictitious "particle" of type *bptype*\ , where *bptype* is either the
   largest atom type in the system, or the type set by the *bptype* flag.
   Any actual existing particles with this atom type will be deleted at
   the beginning of a run. This means you must specify the number of
   types in your system accordingly; usually to be one larger than what
   would normally be the case, e.g. via the :doc:`create\_box <create_box>`
   or by changing the header in your :doc:`data file <read_data>`.  The
   fictitious "bond particles" are inserted at the beginning of the run,
   and serve as placeholders that define the position of the bonds.  This
   allows neighbor lists to be constructed and pairwise interactions to
   be computed in almost the same way as is done for actual particles.
   Because bonds interact only with other bonds, :doc:`pair\_style hybrid <pair_hybrid>` should be used to turn off interactions
   between atom type *bptype* and all other types of atoms.  An error
   will be flagged if :doc:`pair\_style hybrid <pair_hybrid>` is not used.

The optional *exclude* keyword determines if forces are computed
between first neighbor (directly connected) bonds.  For a setting of
*no*\ , first neighbor forces are computed; for *yes* they are not
computed. A setting of *no* cannot be used with the *min* option for
distance calculation because the minimum distance between directly
connected bonds is zero.

Pair style *srp* turns off normalization of thermodynamic properties
by particle number, as if the command :doc:`thermo\_modify norm no <thermo_modify>` had been issued.

The pairwise energy associated with style *srp* is shifted to be zero
at the cutoff distance *rc*\ .


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair styles does not support mixing.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
shift option for the energy of the pair interaction. Note that as
discussed above, the energy term is already shifted to be 0.0 at the
cutoff distance *rc*\ .

The :doc:`pair\_modify <pair_modify>` table option is not relevant for
this pair style.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes global and per-atom information to :doc:`binary restart files <restart>`. Pair srp should be used with :doc:`pair\_style hybrid <pair_hybrid>`, thus the pair\_coeff commands need to be
specified in the input script when reading a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This pair style is part of the USER-MISC package. It is only enabled
if LAMMPS was built with that package. See the Making LAMMPS section
for more info.

This pair style must be used with :doc:`pair\_style hybrid <pair_hybrid>`.

This pair style requires the :doc:`newton <newton>` command to be *on*
for non-bonded interactions.

This pair style is not compatible with :doc:`rigid body integrators <fix_rigid>`

Related commands
""""""""""""""""

:doc:`pair\_style hybrid <pair_hybrid>`, :doc:`pair\_coeff <pair_coeff>`,
:doc:`pair dpd <pair_dpd>`

Default
"""""""

The default keyword value is exclude = yes.


----------


.. _Sirk2:



**(Sirk)** Sirk TW, Sliozberg YR, Brennan JK, Lisal M, Andzelm JW, J
Chem Phys, 136 (13) 134903, 2012.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
