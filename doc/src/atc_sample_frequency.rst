.. index:: fix_modify AtC sample_frequency

fix_modify AtC sample_frequency command
=======================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> sample_frequency <freq>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* sample_frequency = name of the AtC sub-command
* freq = frequency to sample fields in number of steps

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC sample_frequency 10

Description
"""""""""""

Specifies a frequency at which fields are computed for the case where
time filters are being applied.

Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>` and is only relevant when time filters are being used.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None.

