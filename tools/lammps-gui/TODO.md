LAMMPS-GUI TODO list:

# Short term goals

- add CTRL-q hotkey to log windows so you can exit the entire application (add do you really want to? dialog to this)
- add "syntax check" with enabled "-skiprun" flag
- add multi-tab settings dialog where certain properties can be set through customizing the LAMMPS command line
   + select Font
- add list of 5(?) most recently opened/saved files to file dialog (and also write to settings state on exit) (note: must store full path!)

# Long term ideas
- write a "wizard" dialog that can be used for beginners to create an input file template for a few typical use scenarios
- use the "lammps_get_last_thermo" function to get access to thermodynamic data during a run and add plot/graph dialog that can plot one or more of those graphs while the simulation is still running
- possibly also implement a callback interface, so that external programs can be called after thermo data is updated.
- support single stepping, i.e. process input line by line (need to detect continuation chars!) with highlighting active line(s)
- have command text input file in/above status bar where individual commands can be tested. have insert button to copy line into file at the current point
- support text completion as done with lammps-shell
- have context menu for known commands to offer retrieving help by dispatching URL to webbrowser (process index from sphinx for that purpose)
- add a "python" mode, where instead of launching LAMMPS, python is loaded that the LAMMPS python module is made available.
