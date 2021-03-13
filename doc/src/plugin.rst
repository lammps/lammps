.. index:: plugin

plugin command
==============

Syntax
""""""

.. parsed-literal::

   plugin command args

* command = *load* or *unload* or *list*
* args = list of arguments for a particular plugin command

  .. parsed-literal::

     *load* file = load plugin(s) from shared object in *file*
     *unload* style name = unload plugin *name* of style *style*
         *style* = *pair* or *fix* or *command*
     *list* = print a list of currently loaded plugins

Examples
""""""""

.. code-block:: LAMMPS

   plugin load morse2plugin.so
   plugin unload pair morse2/omp
   plugin unload command hello
   plugin list

Description
"""""""""""

The plugin command allows to load (and unload) additional styles and
commands into a LAMMPS binary from so-called dynamic shared object (DSO)
files.  This enables to add new functionality to an existing LAMMPS
binary without having to recompile and link the entire executable.

The *load* command will load and initialize all plugins contained in the
plugin DSO with the given filename.  A message with information the
plugin style and name and more will be printed.  Individual DSO files
may contain multiple plugins.  More details about how to write and
compile the plugin DSO is given in programmer's guide part of the manual
under :doc:`Developer_plugins`.

The *unload* command will remove the given style or the given name from
the list of available styles.  If the plugin style is currently in use,
that style instance will be deleted.

The *list* command will print a list of the loaded plugins and their
styles and names.


Restrictions
""""""""""""

Plugins are currently not available on Windows.

For the loading of plugins to work, the LAMMPS library must be
:ref:`compiled as a shared library <library>`.

Plugins are dependent on the LAMMPS binary interface (ABI)
and particularly the MPI library used. So they are not guaranteed
to work when the plugin was compiled with a different MPI library
or different compilation settings or a different LAMMPS version.
If there is a mismatch the *plugin* command may fail to load the
plugin(s) or data corruption or crashes may happen.


Related commands
""""""""""""""""

none


Default
"""""""

none
