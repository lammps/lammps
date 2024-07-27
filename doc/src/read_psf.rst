.. index:: read_psf

read_psf command
===================

Syntax
""""""

.. code-block:: LAMMPS

   read_psf file

* file = name of psf file to read

Examples
""""""""

.. code-block:: LAMMPS

  read_psf 1hvn-charmm.psf

  write_psf all 1hvn-all.psf

  group solv type OT HT
  group protein_ions subtract all solv
  write_psf protein_ions 1hvn-protein_ions.psf


Description
"""""""""""

Read a file in `Protein Structure Format (psf) <https://www.charmm-gui.org/?doc=lecture&module=pdb&lesson=6>`_. Psf information (segment, residue, atom name, atom type label) is needed by :doc:`write_psf <write_psf>`.

.. admonition:: Segment, Residue, and Atom Name information
   :class: note

   read_psf automatically creates a per-atom custom array (see :doc:`fix_property_atom <fix_property_atom>` and :doc:`compute_property_atom <compute_property_atom>`) called *psf_segment_residue_name* to store SEGMENT, RESIDUE, and NAME information from psf file. If read_psf hasn't been called yet, then this custom array doesn't exist and write_psf will give an error.

.. admonition:: Atom Type Label information
   :class: note

   read_psf automatically creates atom type labels (see :doc:`labelmap <labelmap>`) from psf file. If read_psf hasn't been called yet, then these atom type labels don't exist and write_psf will give an error.


----------

Related commands
""""""""""""""""

:doc:`write_psf <write_psf>`, :doc:`fix_property_atom <fix_property_atom>`, :doc:`compute_property_atom <compute_property_atom>`, :doc:`labelmap <labelmap>`
