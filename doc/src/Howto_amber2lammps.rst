AMBER to LAMMPS Tutorial
========================

**written by Arun Srikanth Sridhar**

.. versionadded:: 11Feb2026

`AMBER2LAMMPS <amber2lammps_>`_ is a Python utility that converts AMBER
topology (``.prmtop``) and coordinate/charge (``.crd``) files into
LAMMPS data and parameter files.  It provides both a command-line
interface and a Python API with built-in validation.  This tutorial
provides a detailed workflow with examples using the `AMBER2LAMMPS
<amber2lammps_>`_ python utility.  The legacy scripts previously
distributed in ``tools/amber2lmp`` have been removed due to their
reliance on Python 2 and lack of maintenance.

This tutorial assumes familiarity with molecular dynamics simulations
and molecular mechanics force fields, including atom typing, bond
topologies, partial charge assignment, and different force field
families (AMBER, CHARMM, GROMOS, OPLS/AA).  For an introduction to force
field concepts, see :doc:`Howto_FFgeneral` before starting this
tutorial.  AMBER offers various force fields for different biomolecules,
including ff19SB/ff14SB (proteins), Bsc1/OL3 (DNA/RNA), GLYCAM06
(carbohydrates), and Lipid17 (lipids), alongside the General AMBER Force
Field (GAFF/GAFF2) for small organic molecules (see
https://ambermd.org/AmberModels.php and Wang, J., et
al., J. Comput. Chem., 25, 1157-1174 (2004) for GAFF force field
development).

.. admonition:: Important Force Field Compatibility Warning:
   :class: warning

   Different force field families (AMBER, CHARMM, GROMOS, OPLS/AA) use
   distinct atom typing schemes, charge assignment methods, and
   parameterization strategies.  **Never mix and match parameters from
   different force field families; this leads to unphysical results and
   simulation failures.** Always use a consistent force field family
   throughout your system.

.. _amber2lammps: https://github.com/askforarun/AMBER2LAMMPS

----------------

.. contents::
   :local:

----------------

What the AMBER2LAMMPS tool does
-------------------------------

This tool helps you run molecular dynamics simulations in **LAMMPS**
when you have your molecular system set up in **AMBER** format.
AMBER2LAMMPS uses AmberTools utilities (antechamber and tleap) to
perform the conversion.

**Typical workflow:**

1. Start with a molecular structure (PDB file or SMILES string to be converted to PDB)
2. Use AMBER tools to create AMBER files (``.prmtop``, ``.crd``)
3. **Convert AMBER files to LAMMPS format** using this tool
4. Run your simulation in LAMMPS

**What gets converted:**

* **``.prmtop``** (topology file): Contains bonds, angles, atom types, and force field parameters -> LAMMPS data and parameter files
* **``.crd``** (coordinates file): Contains atomic positions -> LAMMPS coordinates

**What you get as output:**

* **LAMMPS data file** (e.g., ``data.lammps``): Contains atomic coordinates, box dimensions, and molecular topology
* **LAMMPS parameter file** (e.g., ``parm.lammps``): Contains force field parameters for bonds, angles, and non-bonded interactions

----------------

AMBER2LAMMPS homepage and download
----------------------------------

AMBER2LAMMPS is developed and maintained outside the LAMMPS repository.
Download its source code from its GitHub project page or clone its repository from GitHub:

.. code-block:: bash

   # direct download
   curl -L -o AMBER2LAMMPS-main.tar.gz https://github.com/askforarun/AMBER2LAMMPS/archive/refs/heads/main.tar.gz
   tar -xzvvf AMBER2LAMMPS-main.tar.gz
   cd AMBER2LAMMPS-main

   # clone repository
   git clone https://github.com/askforarun/AMBER2LAMMPS.git
   cd AMBER2LAMMPS

----------------

Requirements
------------

**Platform Compatibility**

AMBER2LAMMPS works on all major platforms:

* **Linux**: Full support (Ubuntu, CentOS, Fedora, etc.)
* **macOS**: Full support (Intel and Apple Silicon)
* **Windows**: Support via WSL2 or Git Bash

**System Requirements**

* **Structure**: PDB file (or a SMILES string that you can convert to PDB).
* **AmberTools utilities**: ``antechamber`` and ``tleap`` to
  build ``.prmtop`` and ``.crd``.
* **Python packages**: ``parmed`` and ``numpy``.
* **LAMMPS**: On your ``PATH`` (``which lmp`` or ``lmp -help``) and built
  with ``MOLECULE``, ``KSPACE``, and ``EXTRA-MOLECULE`` packages.
* **Optional**: Open Babel (``obabel``) if starting from SMILES.

Installing dependencies
-----------------------

AmberTools
^^^^^^^^^^

Install from https://ambermd.org/GetAmber.php#ambertools and activate
the environment:

.. code-block:: bash

   conda activate Ambertools23  # or your AmberTools environment

To see the available force fields, run antechamber -h

Python packages
^^^^^^^^^^^^^^^

Using conda (recommended):

.. code-block:: bash

   conda install -c conda-forge parmed numpy

Using pip:

.. code-block:: bash

   pip install --user parmed numpy

Open Babel (optional, for SMILES to PDB conversion)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using conda (recommended):

.. code-block:: bash

   conda install -c conda-forge openbabel

Using pip:

.. code-block:: bash

   pip install --user openbabel

Using system package manager:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install openbabel

   # Fedora/Red Hat
   sudo dnf install openbabel

   # macOS
   brew install open-babel

----------------

Command Reference
-----------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   * - Argument
     - Required?
     - Description
   * - **Input Files**
     -
     -
   * - ``topology``
     - yes
     - AMBER topology file (``.prmtop``)
   * - ``crd``
     - yes
     - AMBER coordinate file (``.crd``)
   * - **Output Files**
     -
     -
   * - ``data_file``
     - yes
     - Output LAMMPS data filename
   * - ``param_file``
     - yes
     - Output LAMMPS parameter filename
   * - **Options**
     -
     -
   * - ``-b, --buffer``
     - optional
     - Vacuum padding (Angstrom) for the simulation box. Default: ``3.8``.
   * - ``--charge``
     - yes
     - Target net charge (integer). Applies a uniform offset to every
       atom to reach this charge (1e-6 tolerance).
   * -
     -
     - Does not add counterions; add those with the AMBER tools before conversion.
   * - ``--keep-temp``
     - optional
     - Keep temporary files (bonds.txt, angles.txt, dihedrals.txt, pairs.txt) after conversion. Default: ``False``.
   * - ``--verbose``
     - optional
     - Print step-by-step progress, counts, and box size. Default: ``False``.
   * - ``-h, --help``
     - optional
     - Show help message.

Conversion Process
------------------

The converter implements the following sequence (see
``amber_to_lammps.py``):

1. **Input validation**: ``validate_files`` checks that ``topology`` and ``crd`` exist and are readable. The CLI performs this validation automatically, while the Python API requires manual validation if desired.
2. **AMBER topology load**: ParmEd reads atoms, bonds, angles, and
   dihedrals from ``.prmtop``; counts are printed with ``--verbose``.
3. **Atom typing and masses**: Atom types and masses are extracted
   from the topology file.
4. **Coordinates and charges**: Coordinates are read from the
   CRD file and charges from the topology file.
5. **Box creation (``--buffer``)**: A bounding box around the coordinates
   is expanded by the buffer on all sides.
6. **Charge normalization**: Charges are shifted uniformly to achieve the
   target net charge specified by ``--charge``.
7. **Non-bonded parameters**: Lennard-Jones parameters are extracted
   from the topology and become ``pair_coeff`` entries in the parameter file.
8. **Topology terms**: Bonds, angles, and dihedrals are exported via
   ParmEd to temporary files and written to the LAMMPS data/parameter
   files.
9. **Cleanup**: Temporary helper files (``bonds.txt``, ``angles.txt``,
   ``dihedrals.txt``, ``pairs.txt``) are removed unless ``--keep-temp`` is specified.
10. **Verbose diagnostics**: With ``--verbose``, the script reports
    counts, box extents after buffering, and the final total charge.

Charge normalization
^^^^^^^^^^^^^^^^^^^^

**AMBER Charge Schemes:**
AMBER uses various charge methods including:

- **RESP** (Restrained Electrostatic Potential): Derived from quantum
  mechanical electrostatic potential
- **AM1-BCC**: Semi-empirical charges with bond charge corrections
- **CM5** or **CM1A**: Charge models based on atomic charges

In this tutorial, we will use AM1-BCC charges.

**Charge Normalization in AMBER2LAMMPS:**

- Most systems should be overall neutral for PME to converge; small
  residual charges (+/-0.003) from AM1-BCC (or other methods) calculated
  from antechamber are common.
- For neutral molecules, run with ``--charge 0`` so AMBER2LAMMPS applies
  a uniform offset that removes the residual charge; this prevents the
  error from scaling up when the system is replicated in LAMMPS.
- For intentionally charged species (e.g., protonated or deprotonated),
  add counterions in tleap/packmol before conversion. AMBER2LAMMPS never
  adds ions; it only shifts existing charges to your requested total.
- How the flag works:
  - ``--charge 0``: uniform shift makes the summed charge 0 within 1e-6.
  - ``--charge +1`` (or any integer): uniform shift makes the total that integer within 1e-6.
  - The same constant is added to every atom, so relative charge differences are preserved.

**Example:** If your system has net charge +0.003 and you specify
 ``--charge 0``, each atom's charge will be reduced by ``(0.003 /
 number_of_atoms)`` to reach neutrality.

-------------

Workflow Examples
-----------------

Prepare input files (ethanol example)
-------------------------------------

Convert SMILES to PDB (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you start from a SMILES string, generate a PDB with ``obabel``:

.. code-block:: bash

   obabel -:CCO -h -opdb -O ethanol.pdb --gen3d  # adds hydrogens, recommended
   obabel -:c1ccccc1 -h -opdb -O benzene.pdb --gen3d
   obabel -:"CC(=O)OC1=CC=CC=C1C(=O)O" -h -opdb -O aspirin.pdb --gen3d

Generate AMBER topology and coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming ``ethanol.pdb`` is your structure:

1. Generate a MOL2 file with charges:

   .. code-block:: bash

      antechamber -j 4 -at gaff2 -dr yes -fi pdb -fo mol2 \
         -i ethanol.pdb -o ethanol.mol2 -c bcc

   The ``-c bcc`` assigns AM1-BCC charges;
   use the ``-at`` option to choose the force field.

2. Create ``tleap.in``:

   .. code-block:: text

      source leaprc.gaff2
      ETH = loadmol2 ethanol.mol2
      check ETH
      saveamberparm SUS ethanol.prmtop ethanol.crd
      quit

   leaprc.gaff2 loads the gaff2 force field parameters that were used in step 1

3. Run tleap and inspect ``leap.log`` for errors:

   .. code-block:: bash

      tleap -f tleap.in

   Outputs: ``ethanol.prmtop`` and ``ethanol.crd``.


-------------

Basic conversion workflow
-------------------------

LAMMPS input script
^^^^^^^^^^^^^^^^^^^

Save the following LAMMPS commands in a file called ``example_lammps_input.lmp``:

.. code-block:: text

   # LAMMPS Input Script for Converted AMBER System

   units real
   dimension 3
   boundary p p p
   atom_style full

   read_data data.lammps

   pair_style      lj/cut/coul/long 9 9
   bond_style      harmonic
   angle_style     harmonic
   dihedral_style  fourier
   special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333

   include parm.lammps

   thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe
   run 0

The above script includes the parameter file via ``include parm.lammps``
and loads the coordinate data via ``read_data data.lammps``.

CLI usage with LAMMPS execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure amber_to_lammps.py is executable and the input files are accessible.

Run the converter:

.. code-block:: bash

   python3 amber_to_lammps.py data.lammps parm.lammps \
      ethanol.prmtop ethanol.crd --charge 0 --verbose -b 4.5

Outputs: ``data.lammps`` (coordinates/topology) and ``parm.lammps``

Run with LAMMPS:

.. code-block:: bash

   lmp -in example_lammps_input.lmp

Additional CLI examples
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Custom buffer and verbose logging (adds 5 Angstroms of padding)
   python3 amber_to_lammps.py my_data.lammps my_params.lammps \
      ethanol.prmtop ethanol.crd --charge 0 --verbose -b 5.0

   # Custom output names
   python3 amber_to_lammps.py system.data system.parm system.prmtop \
      system.crd --charge 0

   # Minimal output without verbose logging
   python3 amber_to_lammps.py small.data small.parm \
      ethanol.prmtop ethanol.crd --charge 0 -b 3.0

   # Using absolute paths with custom buffer
   python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/topology.prmtop /home/user/amber/coords.crd --charge 0 -b 4.5

   # Keep temporary files for debugging
   python3 amber_to_lammps.py debug_data.lammps debug_parm.lammps molecule.prmtop molecule.crd --charge 0 --keep-temp --verbose

Make sure to rename the files in example_lammps_input.lmp to match the
names of the files generated by the conversion when custom file names
are used.

Python API usage with LAMMPS execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example script uses the :doc:`LAMMPS Python module <Python_head>`
to execute LAMMPS directly from Python.

.. code-block:: python

   from amber_to_lammps import amber2lammps, validate_files
   from lammps import lammps

   validate_files('ethanol.prmtop', 'ethanol.crd')

   amber2lammps(
       data_file='data.lammps',
       param_file='parm.lammps',
       topology='ethanol.prmtop',
       crd='ethanol.crd',
       charge=0,
       buffer=3.8,
       verbose=True,
   )

   lmp = lammps()
   lmp.file('example_lammps_input.lmp')
   lmp.close()

Additional API examples
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Example 1: Basic conversion without validation
   from amber_to_lammps import amber2lammps

   amber2lammps(
       data_file='benzene.data',
       param_file='benzene.parm',
       topology='benzene.prmtop',
       crd='benzene.crd',
       charge=0
   )

   # Example 2: Custom buffer and verbose output
   amber2lammps(
       data_file='ethanol.data',
       param_file='ethanol.parm',
       topology='ethanol.prmtop',
       crd='ethanol.crd',
       charge=0,
       buffer=5.0,
       verbose=True,
       keep_temp=True  # Keep temporary files for inspection
   )

   # Example 3: Batch processing multiple molecules
   from amber_to_lammps import amber2lammps, validate_files

   molecules = ['ethanol', 'benzene', 'aspirin']

   for mol in molecules:
       validate_files(f'{mol}.prmtop', f'{mol}.crd')
       amber2lammps(
           data_file=f'{mol}.data',
           param_file=f'{mol}.parm',
           topology=f'{mol}.prmtop',
           crd=f'{mol}.crd',
           charge=0,
           verbose=True
       )


Validation of AMBER2LAMMPS
--------------------------

AMBER2LAMMPS has been validated against InterMol output. See the project
page for details: https://github.com/askforarun/AMBER2LAMMPS

Getting Help
------------

* **Submit Issues:** https://github.com/askforarun/AMBER2LAMMPS/issues
* **Feature Requests:** Use GitHub Issues or Discussions
* **Questions:** Use GitHub Discussions or Issues

Citation
--------

If you use the AMBER2LAMMPS tool in your research, please cite it as:

**DOI:** `10.5281/zenodo.18114886 <https://doi.org/10.5281/zenodo.18114886>`_

.. raw:: latex

   \clearpage
