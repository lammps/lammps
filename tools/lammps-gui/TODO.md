LAMMPS-GUI TODO list:

# Short term goals

- use qgetenv, qputenv for portability
- rewrite syntax highlighting to be line oriented instead of word oriented.
  handle first part of line based on regular expressions, then advance and only highlight strings and numbers.
  handle "&" continuation and multiline strings with """ like C style comments in Qt docs example
- add CTRL-q hotkey to log windows so you can exit the entire application (add do you really want to? dialog to this)
- add "syntax check" with enabled "-skiprun" flag
- need to handle "label" and "jump" commands from within ?
- switch processing of input to line based commands or?
- switch input file editor to read-only while loop is running
- add Variables menu entry where command line variables can be set

# Long term ideas
- rewrite entire application to either use QtCreator for everything or just build the App and its layout manually
- port to Qt6
- also a rewrite should establish consistent naming conventions. now we have a mix of LAMMPS style, Qt style, and others.
- add option to attach a debugger to the running program (highly non-portable, need customization support in preferences)
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
- support single stepping, i.e. process input line by line (need to detect continuation chars!) with highlighting active line(s)
- have command text input file in/above status bar where individual commands can be tested. have insert button to copy line into file at the current point
- support text completion as done with lammps-shell
- have context menu for known commands to offer retrieving help by dispatching URL to webbrowser (process index from sphinx for that purpose)
- add a "python" mode, where instead of launching LAMMPS, python is loaded that the LAMMPS python module is made available.
