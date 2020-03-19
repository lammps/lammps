.. index:: fix_modify AtC on_the_fly

fix_modify AtC on_the_fly command
=================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> on_the_fly on_the_fly <bond|kernel> <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* on_the_fly = name of the AtC sub-command
* *bond* or *kernel* = specifies on-the-fly calculation of *bond* or *kernel* matrix elements
* *on* or *off* = activate or discontinue on-the-fly mode

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC on_the_fly bond on
   fix_modify AtC on_the_fly kernel
   fix_modify AtC on_the_fly kernel off

Description
"""""""""""

Overrides normal mode of pre-calculating and storing bond pair-to-node a
nd kernel atom-to-node matrices. If activated, it will calculate elements
of these matrices during repeated calls of field computations
(i.e. "on-the-fly") and not store them for future use.  The *on* flag is
optional - if omitted, on_the_fly will be activated for the specified
matrix.  Can be deactivated using the *off* flag.

Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>`.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

By default, on-the-fly calculation is not active (i.e. off). However,
THE code does a memory allocation check to determine if it can store all
needed bond and kernel matrix elements. If this allocation fails,
on-the-fly will be activated.

