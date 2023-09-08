LAMMPS-GUI TODO list:

# Short term goals

- add "syntax check" with enabled "-skiprun" flag
- implement "static" completion for fix/compute/styles/region etc...
- implement "dynamic" completion for variable names, group names, molecule names, compute/dump/fix/region/group IDs
- implement indenting regions for (nested) loops?

# Long term ideas
- rewrite entire application to build the App and its layout manually
- port to Qt6
- also a rewrite should establish consistent naming conventions. now we have a mix of LAMMPS style, Qt style, and others.
- add option to attach a debugger to the running program (highly non-portable, need customization support in preferences)
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
