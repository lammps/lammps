LAMMPS-GUI TODO list:

# Short term goals (v1.x)

- implement indenting regions for (nested) loops?

# Long term ideas (v2.x)
- rewrite entire application to build the App and its layout manually
- port to Qt6 (with compatibility to Qt5?)
- also a rewrite should establish consistent naming conventions. now we have a mix of LAMMPS style, Qt style, and others.
- add option to attach a debugger to the running program (highly non-portable, need customization support in preferences)
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
  (could perhaps use some LLM based KI to look up suggestions for answers?).
