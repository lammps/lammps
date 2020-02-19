.. index:: pair\_style awpmd/cut

pair\_style awpmd/cut command
=============================

Syntax
""""""


.. parsed-literal::

   pair_style awpmd/cut Rc keyword value ...

* Rc = global cutoff, -1 means cutoff of half the shortest box length
* zero or more keyword/value pairs may be appended
* keyword = *hartree* or *dproduct* or *uhf* or *free* or *pbc* or *fix* or *harm* or *ermscale* or *flex\_press*
  
  .. parsed-literal::
  
       *hartree* value = none
       *dproduct* value = none
       *uhf* value = none
       *free* value = none
       *pbc* value = Plen
         Plen = periodic width of electron = -1 or positive value (distance units)
       *fix* value = Flen
         Flen = fixed width of electron = -1 or positive value (distance units)
       *harm* value = width
         width = harmonic width constraint
       *ermscale* value = factor
         factor = scaling between electron mass and width variable mass
       *flex_press* value = none



Examples
""""""""


.. parsed-literal::

   pair_style awpmd/cut -1
   pair_style awpmd/cut 40.0 uhf free
   pair_coeff \* \*
   pair_coeff 2 2 20.0

Description
"""""""""""

This pair style contains an implementation of the Antisymmetrized Wave
Packet Molecular Dynamics (AWPMD) method.  Need citation here.  Need
basic formulas here.  Could be links to other documents.

Rc is the cutoff.

The pair\_style command allows for several optional keywords
to be specified.

The *hartree*\ , *dproduct*\ , and *uhf* keywords specify the form of the
initial trial wave function for the system.  If the *hartree* keyword
is used, then a Hartree multielectron trial wave function is used.  If
the *dproduct* keyword is used, then a trial function which is a
product of two determinants for each spin type is used.  If the *uhf*
keyword is used, then an unrestricted Hartree-Fock trial wave function
is used.

The *free*\ , *pbc*\ , and *fix* keywords specify a width constraint on
the electron wave packets.  If the *free* keyword is specified, then there is no
constraint.  If the *pbc* keyword is used and *Plen* is specified as
-1, then the maximum width is half the shortest box length.  If *Plen*
is a positive value, then the value is the maximum width.  If the
*fix* keyword is used and *Flen* is specified as -1, then electrons
have a constant width that is read from the data file.  If *Flen* is a
positive value, then the constant width for all electrons is set to
*Flen*\ .

The *harm* keyword allow oscillations in the width of the
electron wave packets.  More details are needed.

The *ermscale* keyword specifies a unitless scaling factor
between the electron masses and the width variable mass.  More
details needed.

If the *flex\_press* keyword is used, then a contribution from the
electrons is added to the total virial and pressure of the system.

This potential is designed to be used with :doc:`atom_style wavepacket <atom_style>` definitions, in order to handle the
description of systems with interacting nuclei and explicit electrons.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutoff (distance units)

For *awpmd/cut*\ , the cutoff coefficient is optional.  If it is not
used (as in some of the examples above), the default global value
specified in the pair\_style command is used.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for this pair style.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

These are the defaults for the pair\_style keywords: *hartree* for the
initial wave function, *free* for the wave packet width.
