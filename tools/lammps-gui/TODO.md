LAMMPS-GUI TODO list:

# Short term goals

- add "Help" entry to menu bar. Should open a popup window with a one page description of how to use it. Use HTML or Markdown text.
- display current working directory
- add CTRL-q hotkey to log windows so you can exit the entire application (add do you really want to? dialog to this)
- add "syntax check" with enabled "-skiprun" flag
- add settings dialog where certain properties can be set through customizing the LAMMPS command line
   + enable/disable OpenMP (via suffix), OPT package, select number of OpenMP threads
   + toggle whether captured screen output should include input file echo
   + select Font
   + snapshot image options
- store settings/preferences when changed and read at startup
- add list of 5(?) most recently opened/saved files to file dialog (and also write to settings state on exit) (note: must store full path!)

# Long term ideas
- support single stepping, i.e. process input line by line (need to detect continuation chars!) with highlighting active line(s)
- have command text input file in/above status bar where individual commands can be tested. have insert button to copy line into file at the current point
- support text completion as done with lammps-shell
- have context menu for known commands to offer retrieving help by dispatching URL to webbrowser (process index from sphinx for that purpose)
