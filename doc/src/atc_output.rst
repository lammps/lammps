.. index:: fix_modify AtC output

fix_modify AtC output command
=============================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> output <filename_prefix> <frequency> [text|full_text|binary|vector_components|tensor_components]
   fix_modify <AtC fixID> output index [step|time]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* *output* or *output index* = name of the AtC sub-command
* filename_prefix = prefix for data files (for *output*)
* frequency = frequency of output in time-steps (for *output*)
* optional keywords for *output*:
  - text = creates text output of index, step and nodal variable values for unique nodes
  - full_text = creates text output index, nodal id, step, nodal coordinates and nodal variable values for unique and image nodes
  - binary = creates binary EnSight output
  - vector_components = outputs vectors as scalar components
  - tensor_components = outputs tensor as scalar components (for use with ParaView)
* *step* or *time* = index output by step or by time (for *output index*)


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC output heatFE 100
   fix_modify AtC output hardyFE 1 text tensor_components
   fix_modify AtC output hardyFE 10 text binary tensor_components
   fix_modify AtC output index step


Description
"""""""""""

Creates text and/or binary (EnSight, "gold" format) output of nodal/mesh
data which is transfer/physics specific. Output indexing by step or time
is possible.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix atc command <fix_atc>`

Default
"""""""

No default format. Output indexed by time.
