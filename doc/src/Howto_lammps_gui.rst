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
command to the LAMMPS instance, read the resulting file and show it in
an ``Image Viewer``.  Window.


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
