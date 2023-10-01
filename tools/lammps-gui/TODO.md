LAMMPS-GUI TODO list:

# Short term goals (v1.x)

- implement indenting regions for (nested) loops?
- implement data file manager GUI with the following features:
   - import coordinates and topology via VMD molfile plugins
   - import coordinates and topology from intermol
   - import coordinates and topology from OpenBabel
   - store data internally in a generalized YAML format
   - add/remove columns to per-atom data
   - change atom style for export to data file
   - merge one system to another
   - edit mapping between numeric and symbolic types. create labelmaps.
   - import/export LAMMPS data and molecule files
   - store coordinates internally as unwrapped coordinates
   - recenter coordinates
   - edit box boundaries
   - readjust box to extent of atoms (with or without estimated radius)
   - call to LAMMPS to create geometries from lattices (with/without molecule files) and STL files
   - call to LAMMPS to generate visualizations of geometries
   - edit force field parameters, e.g. apply charmm

# Long term ideas (v2.x)
- rewrite entire application to build the App and its layout manually
- port to Qt6 (with compatibility to Qt5?)
- also a rewrite should establish consistent naming conventions. now we have a mix of LAMMPS style, Qt style, and others.
- add option to attach a debugger to the running program (highly non-portable, need customization support in preferences)
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
  (could perhaps use some LLM based KI to look up suggestions for answers?).
