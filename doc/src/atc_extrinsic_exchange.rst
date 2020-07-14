.. index:: fix_modify AtC extrinsic exchange

fix_modify AtC extrinsic exchange command
=========================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> extrinsic exchange <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* extrinsic exchange = name of the AtC sub-command
* *on* or *off* = set state of energy exchange


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC extrinsic exchange on

Description
"""""""""""

Switches energy exchange between the MD system and the electron system
on or off

Restrictions
""""""""""""

For use only with the two_temperature type of the AtC fix (see
:doc:`fix atc <fix_atc>` command)

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

*on*
