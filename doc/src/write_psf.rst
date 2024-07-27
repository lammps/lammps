.. index:: write_psf

write_psf command
===================

Syntax
""""""

.. code-block:: LAMMPS

   write_psf group-ID file

* group-ID = ID of the group of atoms to be written
* file = name of data file to write psf to

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

Write a file in `Protein Structure Format (PSF) <https://www.charmm-gui.org/?doc=lecture&module=pdb&lesson=6>`_ of the current state of the
simulation. PSF information must have been read in previously by :doc:`read_psf <read_psf>`.

This PSF file is needed together with DCD trajectories created by :doc:`dump dcd <dump>`  for import into visualization software like `ChimeraX <https://www.cgl.ucsf.edu/chimerax/>`_ and `Molecular Nodes for Blender <https://bradyajohnston.github.io/MolecularNodes/>`_.

.. admonition:: Segment, Residue, and Atom Name information
   :class: note

   read_psf automatically creates a per-atom custom array (see :doc:`fix_property_atom <fix_property_atom>` and :doc:`compute_property_atom <compute_property_atom>`) called *psf_segment_residue_name* to store SEGMENT, RESIDUE, and NAME information from PSF file. If read_psf hasn't been called yet, then this custom array doesn't exist and write_psf will give an error.

.. admonition:: Atom Type Label information
   :class: note

   read_psf automatically creates atom type labels (see :doc:`labelmap <labelmap>`) from PSF file. If read_psf hasn't been called yet, then these atom type labels don't exist and write_psf will give an error.

----------

Restrictions
""""""""""""

This command requires inter-processor communication to migrate atoms
before the data file is written.  This means that your system must be
ready to perform a simulation before using this command (force fields
setup, atom masses initialized, etc).

Related commands
""""""""""""""""

:doc:`read_psf <read_psf>`, :doc:`fix_property_atom <fix_property_atom>`, :doc:`compute_property_atom <compute_property_atom>`, :doc:`labelmap <labelmap>`
