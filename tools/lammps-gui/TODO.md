LAMMPS-GUI TODO list:

# Short term goals

- add indicator for when the file in editor is modified (-> status bar?)
- add indicator icon (red/green "LED/lamp") to status bar indicating presense of a "runnable" LAMMPS instance
- add "Help" entry to menu bar. Should open a popup window with a one page description of how to use it. Use HTML or Markdown text.
- add "View" entry to menu bar, where dialog windows may be enabled/disabled (e.g. Render)
- add dialog when exiting asking if file should be saved when it is modified
- add CTRL-q hotkey to log windows so you can exit the entire application (add do you really want to? dialog to this)
- add "render" dialog where a "write_dump image" can be triggered. dialog should offer options for size, zoom, rotation, colors(?)
- add "syntax check" with enabled "-skiprun" flag.
- add settings dialog where certain properties can be set through customizing the LAMMPS command line
   + enable/disable OpenMP (via suffix), select number of OpenMP threads
   + toggle whether captured screen output should include input file echo
   + select Font
- store settings/preferences when changed and read at startup
- add list of recently used files to file dialog (and also write to state on exit) (note: must store full path!)

# Long term ideas
- support single stepping, i.e. process input line by line (need to detect continuation chars!) with highlighting active line(s)
- have command text input file in/above status bar where individual commands can be tested. have insert button to copy line into file at the current point
- support text completion as done with lammps-shell
- have context menu for known commands to offer retrieving help by dispatching URL to webbrowser (process index from sphinx for that purpose)
  (embed html viewer to render pages locally, mathjax??)
