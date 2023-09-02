Using the LAMMPS GUI
====================

This document describes **LAMMPS GUI version 1.5**.

-----

Pre-compiled, ready-to-use LAMMPS GUI executables for Linux (Ubuntu
20.04LTS and compatible or later), macOS (version 11 or later), and
Windows (version 10 or later) :ref:`are available <lammps-gui-install>`.
The source code of the LAMMPS GUI is included in the
``tools/lammps-gui`` folder of the LAMMPS distribution and it can be
compiled alongside LAMMPS.

-----

LAMMPS GUI is a simple graphical text editor that is linked to the
:ref:`LAMMPS C-library interface <lammps_c_api>` and thus can run LAMMPS
directly using the contents of the editor's text buffer as input.  It
can retrieve and display information from LAMMPS while it is running and
is in several ways adapted specifically for editing LAMMPS input files.

This is similar to what people traditionally would do to run LAMMPS
using a command line window: using a regular text editor to edit the
input, run LAMMPS on the input with selected command line flags, and
then extract data from the created files and view them.  That procedure
is quite effective and often required when running LAMMPS on
high-performance computing facilities, or for people proficient in using
the command line, as that would allow them to select tools for the
individual steps they are more comfortable with.  The main benefit of a
GUI application is that it integrates well with graphical desktop
environments and many basic tasks can be done directly from the GUI
without switching to a text console or requiring external programs, let
alone scripts to extract data from the generated output.  This makes it
easier for beginners to get started running simple LAMMPS simulations
and thus very suitable for tutorials on LAMMPS since you only need to
use an learn how to use a single program.  It is designed, however, to
keep the barrier low to switch to a full featured, standalone
programming editor and more sophisticated visualization and analysis
tools later, when running LAMMPS from a command line.

The following text provides a detailed tour of the features and
functionality of the LAMMPS GUI.

-----

Main window
-----------

When LAMMPS GUI starts, it will show the main window with either an
empty buffer or the contents of a file loaded. In the latter case it
may look like the following:

.. image:: JPG/lammps-gui-main.png
   :align: center
   :scale: 50%

There is the menu bar at the top, then the main editor buffer, and a
status bar at the bottom.  The input file contents are shown with line
numbers on the left and the input is colored according to the LAMMPS
input file syntax.  The status bar shows the status of LAMMPS execution
on the left (e.g. "Ready." when idle) and the current working directory
on the right.  The name of the current file in the buffer is shown in
the window title and the text `*modified*` is added in case the buffer
has modifications that are not yet saved to a file.  The size of the
main window will be stored when exiting and restored when starting
again.

Opening Files
^^^^^^^^^^^^^

The LAMMPS GUI application will try to open the first command line
argument as input file, further arguments are ignored.  When no argument
is given, LAMMPS GUI will start with an empty buffer.  Files can also be
opened via the ``File`` menu or by drag-and-drop of a file from a
graphical file manager to the editor window.  Only one file can be open
at a time, so opening a new file with a filled buffer will close this
buffer and - in case the buffer has unsaved modifications - will ask to
either cancel the load, discard the changes, or save them to the file.

Running LAMMPS
^^^^^^^^^^^^^^

From within the LAMMPS GUI main window LAMMPS can be started either from
the ``Run`` menu using the ``Run LAMMPS from Editor Buffer`` entry, by
the hotkey `Ctrl-Enter` (`Command-Enter` on macOS), or by clicking on
the green "Run" button in the status bar.

LAMMPS runs in a separate thread, so the GUI stays responsive and is is
able to interact with the running calculation and access its data.  It
is important to note, that running LAMMPS this way is using the contents
of the input buffer for the run (via the
:cpp:func:`lammps_commands_string()` function of the LAMMPS C-library
interface), and **not** the file it was read from.  Thus, if there are
unsaved changes in the buffer, they *will* be used.  As an alternative,
it is also possible to start LAMMPS by reading the contents of the file
from the ``Run LAMMPS from File`` menu entry or with `Ctrl-Shift-Enter`.
This option may be required in some rare cases where the input uses some
functionality that is not compatible with running LAMMPS from a string
buffer.  For consistency, any unsaved changes in the buffer must be
either saved to the file or undone before LAMMPS can be run from a file.

.. image:: JPG/lammps-gui-running.png
   :align: center
   :scale: 75%

While LAMMPS is running, the contents of the status bar change: on the
left side there is a text indicating that LAMMPS is running, which will
contain the selected number of threads, if thread-parallel acceleration
was selected in the ``Preferences`` dialog.  On the right side, a
progress bar is shown that displays the estimated progress on the
current :doc:`run command <run>`.

Also, the line number of the currently executed command will be
highlighted in green.

.. image:: JPG/lammps-gui-run-highlight.png
   :align: center
   :scale: 75%


In case of an error (in the example below the command :doc:`label
<label>` was incorrectly capitalized as "Label"), an error message
dialog will be shown and the line of the input where the error was
triggered will be highlighted.  The state of LAMMPS as shown in the
status bar will be set to "Failed." instead of "Ready."

.. image:: JPG/lammps-gui-run-error.png
   :align: center
   :scale: 75%

Additionally,  up to three windows will open during a run:

- a log window with the captured screen output
- a chart window with a line graph created from the thermodynamic output of the run
- a slide show window with images created by a :doc:`dump image command <dump_image>`

More information on those windows and how to adjust their behavior and
contents is below.

An active LAMMPS run can be stopped cleanly by using either the ``Stop
LAMMPS`` entry in the ``Run`` menu, the hotkey `Ctrl-/` (`Command-/` on
macOS), or by clicking on the red button in the status bar.  This will
cause that the running LAMMPS process to complete the current iteration
and then complete the processing the input while skipping all run or
minimize commands.  This is equivalent to the command :doc:`timer
timeout 0 <timer>` and implemented by calling the
:cpp:func:`lammps_force_timeout()` function of the LAMMPS C-library
interface.  Please see the corresponding documentation pages to
understand the implications of this feature.

Log Window
----------

By default, when starting a run, a "Log Window" will open that displays
the current screen output of the LAMMPS calculation, that would normally
be seen in the command line window, as shown below.

.. image:: JPG/lammps-gui-log.png
   :align: center
   :scale: 50%

LAMMPS GUI captures the screen output and updates the log window
regularly during a run with it as it is generated.

By default, there will be a new window for each run, so that it is
possible to visually compare outputs from different runs.  It is also
possible to change the behavior of LAMMPS GUI in the preferences dialog
to *replace* an existing log window for a new run or to not show the log
window by default.  It is also possible to show or hide the current log
window from the ``View`` menu.

The text in the log window is read-only and cannot be modified, but
editor commands to select and copy all or parts of the text can be used.
The "Select All" and "Copy" functions are also available via a context
menu by clicking with the right mouse button.

Chart Window
------------

By default, when starting a run, a "Chart Window" will open that
displays a line chart of thermodynamic output of the LAMMPS calculation
as shown below.

.. image:: JPG/lammps-gui-chart.png
   :align: center
   :scale: 50%

The drop down menu on the top right allows to select between the
different properties that are computed and written to the output.  Only
one property can be shown at a time.  These charts will be updated with
new data as the run progresses, so they can be used to visually monitor
the evolution of the available properties.  From the ``File`` menu on
the top left, it is possible to save an image of the currently displayed
chart or export the data in either plain text columns (as usable for
plotting tools like `gnuplot <http://www.gnuplot.info/>`_ or `grace
<https://plasma-gate.weizmann.ac.il/Grace/>`_), or as CSV data which can
be imported for further processing with Microsoft Excel or `pandas
<https://pandas.pydata.org/>`_

Data from multiple successive run commands will be combined into a
single data set unless the format, number, or names of output columns
are changed with a :doc:`thermo_style <thermo_style>` or
:doc:`thermo_modify <thermo_modify>` or the current time step is reset
with :doc:`reset_timestep <reset_timestep>` or if a :doc:`clear <clear>`
command is issued.

Image Slide Show
----------------

By default, in case the LAMMPS input contains a :doc:`dump image
<dump_image>` command, a "Slide Show" window will open which loads and
displays the images created by LAMMPS as they are written.

.. image:: JPG/lammps-gui-slideshow.png
   :align: center
   :scale: 50%

The various buttons at the bottom right of the window allow to either
single step through the list of images or play an animation (as a
continuous loop or once from first to last).  It is also possible to
zoom in or zoom out if the displayed image.

Variable Info
-------------

During a run, it may be of interest to monitor the value of variables,
for example to monitor the progress of loops.  This can be done via
enabling the "Variables Window" in the ``View`` menu or by using the
`Ctrl-Shift-W` hotkey.  This will show info similar to the :doc:`info
variables <info>` command in a separate window as shown below.

.. image:: JPG/lammps-gui-variable-info.png
   :align: center
   :scale: 75%

Like the log and chart windows, its content is continuously updated
during a run, and will show "(none)" if there are no variables defined.
Please note that it is also possible to *set* :doc:`index style
variables <variable>`, that would normally be set via command line flags,
via the "Set Variables..." dialog from the ``Run`` menu.

Viewing Snapshot Images
-----------------------

By selecting the ``Create Image`` entry in the ``Run`` menu, by hitting
the `Ctrl-I` (`Command-I` on macOS) hotkey, or by clicking on the
"palette" button in the status bar, LAMMPS GUI will send a custom
:doc:`write_dump image <dump_image>` command and read the resulting
snapshot image with the current state of the system into an image viewer
window.  This functionality is not available *during* an ongoing run.

When possible, LAMMPS GUI will try to detect which elements the atoms
correspond to (via their mass) and then colorize them in the image
accordingly.  Otherwise the default predefined sequence of colors is
assigned to the different atom types.

.. image:: JPG/lammps-gui-image.png
   :align: center
   :scale: 50%

The default image size, some default image quality settings, the view
style and some colors can be changed in the ``Preferences`` dialog
window.  From the image viewer window further adjustments can be made:
actual image size, high-quality rendering, anti-aliasing, view style,
display of box or axes, zoom factor.  The view on the system can be
rotated horizontally and vertically, and it is possible to only display
the atoms within a group defined in the input (default is "all").  After
each change, the image is rendered again and the display updated.  The
small palette icon on the top left will be colored while LAMMPS is
running to render the new image and it will be grayed out again, when it
is done.  When there are many items to show and high quality images with
anti-aliasing are requested, re-rendering can take several seconds.
From the ``File`` menu of the image window, the shown image can be saved
to a file or copied into the cut-n-paste buffer for pasting into another
application.


Editor Functions
----------------

The editor has most the usual functionality that similar programs have:
text selection via mouse or with cursor moves while holding the Shift
key, Cut (`Ctrl-X`), Copy (`Ctrl-C`), Paste (`Ctrl-V`), Undo (`Ctrl-Z`),
Redo (`Ctrl-Shift-Z`), Select All (`Ctrl-A`).  All of these editing
functions are available via the indicated hotkeys.  When trying to exit
the editor with a modified buffer, a dialog will pop up asking whether
to cancel the quit, or don't save or save the buffer's contents to a
file.

Context Specific Word Completion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, LAMMPS GUI will display a popup window with possible
completions for commands or styles after 3 characters of a word have
been typed. The word can then be completed through selecting an entry by
scrolling down with the cursor keys and selecting with the 'Enter' key
or by clicking on the entry with the mouse.  The automatic completion
popup can be disabled in the ``Preferences`` dialog, but the completion
can still be requested manually by either hitting 'Shift-TAB' key or by
right-clicking with the mouse and selecting the option from the context
menu.  Most of the completion information is taken from the LAMMPS
instance and thus it will be adjusted to only show options available
that have been enabled while compiling LAMMPS, however that excludes
accelerated styles and commands.  Only non-suffix versions are shown.

Line Reformatting
^^^^^^^^^^^^^^^^^

The editor supports reformatting lines according to the syntax in order
to have consistently aligned lines.  This primarily means to add padding
to commands, type specifiers, IDs and names.  This reformatting is
performed by default when hitting the 'Enter' key to start a new line.
This feature can be turned off in the ``Preferences`` dialog, but it can
still be manually performed by hitting the 'TAB' key.

Internally this functionality is achieved by splitting the line into
"words" and then putting it back together with padding added where
the context can be detected; otherwise a single blank is used.

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
letter, that is `Alt-F` will activate the ``File`` menu.  For the corresponding
activated sub-menus, also the underlined letter, together with the `Alt` key can
be used to select entries instead of the using mouse.

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
  will appear to either cancel the quit, to save, or to not save the
  edited file.

In addition, up to 5 recent file names will be listed after the ``Open``
entry that allows to re-open those recent files. This list is stored when
quitting and recovered when starting again.

Edit
^^^^

The ``Edit`` menu offers the usual editor functions like ``Undo``,
``Redo``, ``Cut``, ``Copy``, ``Paste``, but also offers to open the
``Preferences`` dialog (hotkey `Ctrl-P`) and to delete all stored
preferences so they will be reset to their default values.

Run
^^^

The ``Run`` menu allows to start and stop a LAMMPS process.  Rather than
calling the LAMMPS executable as a separate executable, the LAMMPS GUI
is linked to the LAMMPS library and thus can run LAMMPS internally
through the :ref:`LAMMPS C-library interface <lammps_c_api>`.

Specifically, a LAMMPS instance will be created by calling
:cpp:func:`lammps_open_no_mpi` and then the buffer contents are run by
calling :cpp:func:`lammps_commands_string`.  Certain commands and
features are only available, after a LAMMPS instance is created.  Its
presence is indicated by a small LAMMPS ``L`` logo in the status bar at
the bottom left of the main window.  As an alternative, it is also
possible to run LAMMPS using the contents of the edited file by reading
the file.  This is mainly provided as a fallback option in case the
input uses some feature that is not available when running from a string
buffer.

The LAMMPS calculation will be run in a concurrent thread so that the
GUI will stay responsive and will be updated during the run.  This can
be used to tell the running LAMMPS instance to stop at the next
timestep.  The ``Stop LAMMPS`` entry will do this by calling
:cpp:func:`lammps_force_timeout`, which is equivalent to a :doc:`timer
timeout 0 <timer>` command.

The ``Set Variables...`` entry will open a dialog box where :doc:`index
style variables <variable>` can be set. Those variables will be passed
to the LAMMPS instance when it is created and are thus set *before* a
run is started.

.. image:: JPG/lammps-gui-variables.png
   :align: center
   :scale: 75%

The ``Set Variables`` dialog will be pre-populated with entries that are
set as index variables in the input and any variables that are used but
not defined as far as the built-in parser can detect them.  New rows for
additional variables can be added through the ``Add Row`` button and
existing rows may be deleted by clicking on the ``X`` icons on the right.

The ``Create Image`` entry will send a :doc:`dump image <dump_image>`
command to the LAMMPS instance, read the resulting file, and show it in
an ``Image Viewer`` window.

The ``View in OVITO`` entry will launch `OVITO <https://ovito.org>`_
with a :doc:`data file <write_data>` of the current state of the system.
This option is only available, if the LAMMPS GUI can find the OVITO
executable in the system path.

The ``View in VMD`` entry will instead launch VMD, also to load a
:doc:`data file <write_data>` of the current state of the system.  This
option is only available, if the LAMMPS GUI can find the VMD executable
in the system path.

View
^^^^

The ``View`` menu offers to show or hide the additional windows with log
output, charts, slide show, variables, or snapshot images.  The default
settings for those can be changed in the ``Preferences dialog``.

About
^^^^^

The ``About`` menu finally offers a couple of dialog windows and an
option to launch the LAMMPS online documentation in a web browser.  The
``About LAMMPS`` entry displays a dialog with a summary of the
configuration settings of the LAMMPS library in use and the version
number of LAMMPS GUI itself.  The ``Quick Help`` displays a dialog with
a minimal description of LAMMPS GUI.  The ``LAMMPS GUI Howto`` entry
will open this documentation page from the online documentation in a web
browser window.  And ``LAMMPS Manual`` will open the main page of the
LAMMPS documentation in the web browser.

-----

Preferences
-----------

The ``Preferences`` dialog allows to customize some of the behavior
and looks of the LAMMPS GUI application.  The settings are grouped
and each group is displayed within a tab.

.. |guiprefs1| image:: JPG/lammps-gui-prefs-general.png
   :width: 24%

.. |guiprefs2| image:: JPG/lammps-gui-prefs-accel.png
   :width: 24%

.. |guiprefs3| image:: JPG/lammps-gui-prefs-image.png
   :width: 24%

.. |guiprefs4| image:: JPG/lammps-gui-prefs-editor.png
   :width: 24%

|guiprefs1|  |guiprefs2|  |guiprefs3|  |guiprefs4|

General Settings:
^^^^^^^^^^^^^^^^^

- *Echo input to log:* when checked, all input commands, including
  variable expansions, will be echoed to the log window. This is
  equivalent to using `-echo screen` at the command line.  There is no
  log *file* produced by default, since LAMMPS GUI uses `-log none`.
- *Include citation details:* when checked full citation info will be
  included to the log window.  This is equivalent to using `-cite
  screen` on the command line.
- *Show log window by default:* when checked, the screen output of a
  LAMMPS run will be collected in a log window during the run
- *Show chart window by default:* when checked, the thermodynamic
  output of a LAMMPS run will be collected and displayed in a chart
  window as line graphs.
- *Show slide show window by default:* when checked, a slide show
  window will be shown with images from a dump image command, if
  present, in the LAMMPS input.
- *Replace log window on new run:* when checked, an existing log
  window will be replaced on a new LAMMPS run, otherwise each run will
  create a new log window.
- *Replace chart window on new run:* when checked, an existing chart
  window will be replaced on a new LAMMPS run, otherwise each run will
  create a new chart window.
- *Replace image window on new render:* when checked, an existing
  chart window will be replaced when a new snapshot image is requested,
  otherwise each command will create a new image window.
- *Path to LAMMPS Shared Library File:* this options is only available
  when LAMMPS GUI was compiled to load the LAMMPS library at run time
  instead of being linked to it directly.  With the ``Browse..`` button
  or by changing the text, a different shared library file with a
  different compilation of LAMMPS with different settings or from a
  different version can be loaded.  After this setting was changed,
  LAMMPS GUI needs to be re-launched.
- *Select Default Font:* Opens a font selection dialog where the type
  and size for the default font (used for everything but the editor and
  log) of the application can be set.
- *Select Text Font:* Opens a font selection dialog where the type and
  size for the text editor and log font of the application can be set.

Accelerators:
^^^^^^^^^^^^^

This tab enables to select which accelerator package is used and is
equivalent to using the `-suffix` and `-package` flags on the command
line.  Only settings supported by the LAMMPS library and local hardware
are available.  The `Number of threads` field allows to set the maximum
number of threads for the accelerator packages that use threads.

Snapshot Image:
^^^^^^^^^^^^^^^

This tab allows to set some defaults for the snapshot images displayed
in the ``Image Viewer`` window, like its dimensions and the zoom factor
applied.  The *Antialias* switch requests to render images with twice
the number of pixels for width and height and then smoothly scales the
image back to the requested size.  This produces higher quality images
with smoother edges at the expense of requiring more CPU time to render
the image.  The *HQ Image mode* option turns on using a screen space
ambient occlusion mode (SSAO) when rendering images.  This is also more
time consuming, but produces a more 'spatial' representation of the
system.  The *VDW Style* checkbox selects whether atoms are represented
by space filling spheres when checked or by smaller spheres and stick.
Finally there are a couple of drop down lists to select the background
and box color.

Editor Settings:
^^^^^^^^^^^^^^^^

This tab allows to tweak some settings of the editor window.  Specifically
the amount of padding to be added to LAMMPS commands, types or type ranges,
IDs (e.g. for fixes), and names (e.g. for groups).  The value set is the
minimum width for the text element and it can be chosen in the range between
1 and 32.

The following two settings allow to enable or disable the automatic
reformatting on hitting the 'Enter' key and the automatic display of the
completion popup window.

-----------

Hotkeys
-------

Almost all functionality is accessible from the menu or via hotkeys.
The following hotkeys are available (On macOS use the Command key
instead of Ctrl/Control).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Hotkey
     - Function
     - Hotkey
     - Function
     - Hotkey
     - Function
   * - Ctrl+N
     - New File
     - Ctrl+Z
     - Undo edit
     - Ctrl+Enter
     - Run Input
   * - Ctrl+O
     - Open File
     - Ctrl+Shift+Z
     - Redo edit
     - Ctrl+/
     - Stop Active Run
   * - Ctrl+S
     - Save File
     - Ctrl+C
     - Copy text
     - Ctrl+Shift+V
     - Set Variables
   * - Ctrl+Shift+S
     - Save File As
     - Ctrl+X
     - Cut text
     - Ctrl+I
     - Snapshot Image
   * - Ctrl+Q
     - Quit
     - Ctrl+V
     - Paste text
     - Ctrl+L
     - Slide Show
   * - Ctrl+W
     - Close Window
     - Ctrl+A
     - Select All
     - Ctrl+P
     - Preferences
   * - Ctrl+Shift+A
     - About LAMMPS
     - Ctrl+Shift+H
     - Quick Help
     - Ctrl+Shift+G
     - LAMMPS GUI Howto
   * - Ctrl+Shift+M
     - LAMMPS Manual
     - Ctrl+?
     - Context Help
     - Ctrl+Shift+W
     - Show Variables
   * - Ctrl+Shift+Enter
     - Run File
     - TAB
     - Reformat line
     - Shift+TAB
     - Show Completions

Further editing keybindings `are documented with the Qt documentation
<https://doc.qt.io/qt-5/qplaintextedit.html#editing-key-bindings>`_.  In
case of conflicts the list above takes precedence.
