Using the LAMMPS GUI
====================

LAMMPS GUI is essentially a simple graphical text editor that is linked
to the :ref:`LAMMPS C-library interface <lammps_c_api>` and thus can run
LAMMPS directly using the contents of the editor's text buffer as input.

This is similar to what people usually would do to run LAMMPS using a
regular text editor to edit the input and run the necessary command,
possibly including the text editor, from a command line terminal window.
That is quite effective when running LAMMPS on high-performance
computing facilities and when you are very proficient with using the
command line.  The main benefit of a GUI application is that this
integrates very well with graphical desktop environments and many basic
tasks can be done directly from within the GUI without switching to a
text console and requiring external programs or scripts to extract data
from the generated output.  This makes it easier for beginners to get
started running simple LAMMPS simulation and thus very suitable for
tutorials on LAMMPS and then makes it easier to switch to a full
featured text editor and more sophisticated visualization and analysis
tools.

The following text provides a detailed tour of the features and
functionality of the LAMMPS GUI.

Main window
-----------

When LAMMPS GUI starts, it will show the main window with either an
empty buffer, or have a file loaded. In the latter case it may look like
the following:

.. image:: JPG/lammps-gui-main.png
   :align: center
   :scale: 50%

There is the menu bar at the top, then the main editor buffer with the
input file contents in the center with line numbers on the left and the
input colored according to some LAMMPS syntax rules.  At the bottom is
the status bar, which shows the status of LAMMPS execution on the right
("Ready." when idle) and the current working directory on the left.
The size of the main window will be stored when exiting and restored
when starting again.  The name of the current file in the buffer is
shown in the window title and the text `*modified*` is added in case
the buffer has modifications that are not yet saved to a file.

Editor Functions
^^^^^^^^^^^^^^^^

The editor has the usual functionality that similar programs have: text
selection via mouse or with cursor moves while holding the Shift key,
Cut, Copy, Paste, Undo, Redo.  All of these editing functions are available
via hotkeys.  When trying to exit the editor with a modified buffer, a
dialog will pop up asking whether to cancel the quit, or don't save or
save the buffer's contents to a file.

Context Specific Help
^^^^^^^^^^^^^^^^^^^^^

.. image:: JPG/lammps-gui-popup-help.png
   :align: center
   :scale: 50%

A unique feature of the LAMMPS GUI is the option to look up the
documentation for the command in the current line.  This can be achieved
by either clicking the right mouse button or by using the `Ctrl-?`
hotkey.  When clicking the mouse there are additional entries in the
context menu that will open the corresponding documentation page in the
online LAMMPS documentation.  When using the hotkey, the first of those
entries will be chosen directly.

Menu
----

The menu bar the entries ``File``, ``Edit``, ``Run``, ``View``, and ``About``.
Instead of using the mouse to click on them, the individual menus can also
be activated by hitting the `Alt` key together with the corresponding underlined
letter, that is `Alt-f` will activate the ``File`` menu.  For the corresponding
activated sub-menus, also the underlined letter, together with the `Alt` key can
be used to select instead of the mouse.

File
^^^^

The ``File`` menu offers the usual options:

- ``New`` will clear the current buffer and reset the file name to ``*unknown*``
- ``Open`` will open a dialog to select a new file
- ``Save`` will save the current file; if the file name is ``*unknown*``
  a dialog will open to select a new file name
- ``Save As`` will open a dialog to select and new file name and save
  the buffer to it
- ``Quit`` will exit LAMMPS GUI. If there are unsaved changes, a dialog
  will appear to either cancel the quit, save or don't save the file.

In addition, up to 5 recent file names will be listed after the ``Open``
entry that allows to re-open recent files. This list is stored when
quitting and recovered when starting again.

Edit
^^^^

The ``Edit`` menu offers the usual editor functions like ``Undo``,
``Redo``, ``Cut``, ``Copy``, ``Paste``, but also offers to open the
``Preferences`` dialog and to delete all stored preferences so they
will be reset to their defaults.

Run
^^^

The ``Run`` menu allows to start and stop a LAMMPS process.  Rather than
calling the LAMMPS executable as a separate executable, the LAMMPS GUI
is linked to the LAMMPS library and thus can run LAMMPS internally
through the :ref:`LAMMPS C-library interface <lammps_c_api>`.
Specifically, a LAMMPS instance will be created by calling
:cpp:func:`lammps_open_no_mpi` and then the buffer contents run by
calling :cpp:func:`lammps_commands_string`.  Certain commands and
features are only available, after a LAMMPS instance is created.  Its
presence is indicated by a small LAMMPS ``L`` logo in the status bar at
the bottom left of the main window.

The LAMMPS calculation will be run in a concurrent thread so that the
GUI will stay responsive and will be updated during the run.  This can
be used to tell the running LAMMPS instance to stop at the next
timestep.  The ``Stop LAMMPS`` entry will do this by calling
:cpp:func:`lammps_force_timeout`, which is equivalent to a :doc:`timer
timeout 0 <timer>` command.

The ``Set Variables`` entry will open a dialog box where :doc:`index style variables <variable>`
can be set. Those variables will be passed to the LAMMPS instance when
it is created and are thus set *before* a run is started.

.. image:: JPG/lammps-gui-variables.png
   :align: center
   :scale: 75%

The ``Set Variables`` dialog will be pre-populated with entries that are
set as index variables in the input and any variables that are used but
not defined as far as the built-in parser can detect them.  New rows for
additional variables can be added through the ``Add Row`` button and
existing rows deleted by clicking on the ``X`` icons on the right.

The ``View Image`` entry will send a :doc:`dump image <dump_image>`
command to the LAMMPS instance, read the resulting file, and show it in
an ``Image Viewer`` window.

View
^^^^

The ``View`` menu offers to show or hide the three optional windows
with log output, graphs, or images.  The default settings for those
can be changed in the ``Preferences dialog``.

About
^^^^^

The ``About`` menu finally offers a couple of dialog windows and an
option to launch the LAMMPS online documentation in a web browser.  The
``About LAMMPS GUI`` entry displays a dialog with a summary of the
configuration settings of the LAMMPS library in use and the version
number of LAMMPS GUI itself.  The ``Quick Help`` displays a dialog with
a minimal description of LAMMPS GUI.  And ``LAMMPS Manual`` will open
the main page of this LAMMPS documentation at https://docs.lammps.org/.

Preferences
-----------

The ``Preferences`` dialog allows to customize some of the behavior
and looks of the LAMMPS GUI application.  The settings are grouped
and each group is displayed within a tab.

.. |guiprefs1| image:: JPG/lammps-gui-prefs-general.png
   :width: 32%

.. |guiprefs2| image:: JPG/lammps-gui-prefs-accel.png
   :width: 32%

.. |guiprefs3| image:: JPG/lammps-gui-prefs-image.png
   :width: 32%

|guiprefs1|  |guiprefs2|  |guiprefs3|

General Settings:
^^^^^^^^^^^^^^^^^

- *Echo input to log:* when checked, all input commands, including
  variable expansions, will be echoed to the log window. This is
  equivalent to using `-echo screen` at the command line.  There is no
  log *file* produced since it always uses `-log none`.
- *Include citation details:* when checked full citation info will be
  included to the log window.  This is equivalent to using `-cite
  screen` on the command line.
- *Show log window by default:* when checked, the screen output of a
  LAMMPS run will be collected in a log window during the run
- *Show chart window by default:* when checked, the thermodynamic
  output of a LAMMPS run will be collected and displayed in a chart
  window as line graphs.
- *Replace log window on new run:* when checked, an existing log
  window will be replaced on a new LAMMPS run, otherwise each run will
  create a new log window.
- *Replace chart window on new run:* when checked, an existing chart
  window will be replaced on a new LAMMPS run, otherwise each run will
  create a new chart window.
- *Replace image window on new render:* when checked, an existing
  chart window will be replaced when a new snapshot image is requested,
  otherwise each command will create a new image window.
- *Select Default Font:* Opens a font selection dialog where the type
  and size for the default font (used for everthing but the editor and
  log) of the application can be set.
- *Select Text Font:* Opens a font selection dialog where the type and
  size for the text editor and log font of the application can be set.

Accelerators:
^^^^^^^^^^^^^

This tab enables to select accelerator settings and is equivalent to
using the `-suffix` and `-package` flags on the command line.  Only
settings supported by the LAMMPS library and local hardware are
available.  The `Number of threads` field allows to set the maximum
number of threads for the accelerator packages that use threads.

Snapshot Image:
^^^^^^^^^^^^^^^

This tab allows to set some defaults for the snapshot images displayed
in the ``Image Viewer`` window, like its dimensions and the zoom factor
applied.  The *Antialias* switch requests to render images at an up to 4
times larger image size and then uses a bicubic scaling algorithm to
rescale them back to the requested size.  This produces higher quality
images with smoother edges at the expense of requiring more CPU time to
render the image.  The *HQ Image mode* option turns on using a screen
space ambient occlusion mode when rendering images.  This is more time
consuming, but produces a more 'spatial' representation of the system.


Hotkeys
-------

Almost all functionality is accessible from the menu or via hotkeys.
The following hotkeys are available (On macOS use the Command key
instead of Ctrl (aka Control)).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Hotkey
     - Function
     - Hotkey
     - Function
     - Hotkey
     - Function
     - Hotkey
     - Function
   * - Ctrl+N
     - New File
     - Ctrl+Z
     - Undo edit
     - Ctrl+V
     - Paste text
     - Ctrl+Q
     - Quit (Main Window only)
   * - Ctrl+O
     - Open File
     - Ctrl+Shift+Z
     - Redo edit
     - Ctrl+Enter
     - Run LAMMPS
     - Ctrl+W
     - Close (Log and Image Window only)
   * - CTRL+S
     - Save File
     - Ctrl+C
     - Copy text
     - Ctrl+/
     - Stop Active Run
     - Ctrl+P
     - Preferences
   * - Ctrl+Shift+S
     - Save File As
     - Ctrl+X
     - Cut text
     - Ctrl+I
     - Create Snapshot Image
     - Ctrl+Shift+/
     - Quick Help

Further editing keybindings `are documented with the Qt documentation
<https://doc.qt.io/qt-5/qplaintextedit.html#editing-key-bindings>`_.  In
case of conflicts the list above takes precedence.
