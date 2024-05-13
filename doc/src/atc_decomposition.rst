.. index:: fix_modify AtC decomposition

fix_modify AtC decomposition command
====================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> decomposition <type>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* decomposition = name of the AtC sub-command
* type =  *replicated_memory* or *distributed_memory*


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC decomposition distributed_memory

Description
"""""""""""

Command for assigning the distribution of work and memory for parallel
runs.  With *replicated_memory* the nodal information is replicated on
each processor, and with *distributed_memory* only the owned nodal
information kept on each processor.  The *replicated_memory* option
is most appropriate for simulations were the number of nodes is much
smaller than the number of atoms.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

replicated_memory
