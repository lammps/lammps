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
