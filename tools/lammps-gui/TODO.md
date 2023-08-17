LAMMPS-GUI TODO list:

# Short term goals

- rewrite syntax highlighting to be line oriented instead of word oriented.
  handle first part of line based on regular expressions, then advance and only highlight strings and numbers.
  handle "&" continuation and multiline strings with """ like C style comments in Qt docs example
- add "syntax check" with enabled "-skiprun" flag
- need to handle "label" and "jump" commands from within ?
- switch input file editor to read-only while loop is running
- triple quoted heredocs don't work with lammps_commands_string()

# Long term ideas
- rewrite entire application to build the App and its layout manually
- port to Qt6
- also a rewrite should establish consistent naming conventions. now we have a mix of LAMMPS style, Qt style, and others.
- add option to attach a debugger to the running program (highly non-portable, need customization support in preferences)
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
- support single stepping, i.e. process input line by line (need to detect continuation chars!) with highlighting active line(s)
- have command text input line in/above status bar where individual commands can be tested. have insert button to copy line into file at the current point
- support text completion as done with lammps-shell
