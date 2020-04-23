.. index:: fix_modify AtC filter type

fix_modify AtC filter type command
===================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> filter type <exponential|step|no_filter>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* filter type = name of the AtC sub-command
* *exponential* or *step* or *no_filter* = select type of filter

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC filter type exponential

Description
"""""""""""

Specifies the type of time filter used.

Restrictions
""""""""""""

Only for use with these specific transfers: thermal, two_temperature

Related AtC commands
""""""""""""""""""""
- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC filter <atc_time_filter>`
- :doc:`fix_modify AtC filter scale <atc_filter_scale>`

Default
"""""""

None.
