Writing plugins
---------------

Plugins provide a mechanism to add functionality to a LAMMPS executable
without recompiling LAMMPS.  The functionality for this and the
:doc:`plugin command <plugin>` are implemented in the
:ref:`PLUGIN package <PKG-PLUGIN>` which must be installed to use plugins.

Plugins use the operating system's capability to load dynamic shared
object (DSO) files in a way similar shared libraries and then reference
specific functions in those DSOs.  Any DSO file with plugins has to include
an initialization function with a specific name, "lammpsplugin_init", that
has to follow specific rules described below.  When loading the DSO with
the "plugin" command, this function is looked up and called and will then
register the contained plugin(s) with LAMMPS.

From the programmer perspective this can work because of the object
oriented design of LAMMPS where all pair style commands are derived from
the class Pair, all fix style commands from the class Fix and so on and
usually only functions present in those base classes are called
directly.  When a :doc:`pair_style` command or :doc:`fix` command is
issued a new instance of such a derived class is created.  This is done
by a so-called factory function which is mapped to the style name.  Thus
when, for example, the LAMMPS processes the command ``pair_style lj/cut
2.5``, LAMMPS will look up the factory function for creating the
``PairLJCut`` class and then execute it.  The return value of that
function is a ``Pair *`` pointer and the pointer will be assigned to the
location for the currently active pair style.

A DSO file with a plugin thus has to implement such a factory function
and register it with LAMMPS so that it gets added to the map of available
styles of the given category.  To register a plugin with LAMMPS an
initialization function has to be present in the DSO file called
``lammpsplugin_init`` which is called with three ``void *`` arguments:
a pointer to the current LAMMPS instance, a pointer to the opened DSO
handle, and a pointer to the registration function.  The registration
function takes two arguments: a pointer to a ``lammpsplugin_t`` struct
with information about the plugin and a pointer to the current LAMMPS
instance.  Please see below for an example of how the registration is
done.

Members of ``lammpsplugin_t``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Member
     - Description
   * - version
     - LAMMPS Version string the plugin was compiled for
   * - style
     - Style of the plugin (pair, bond, fix, command, etc.)
   * - name
     - Name of the plugin style
   * - info
     - String with information about the plugin
   * - author
     - String with the name and email of the author
   * - creator.v1
     - Pointer to factory function for pair, bond, angle, dihedral, improper or command styles
   * - creator.v2
     - Pointer to factory function for compute, fix, or region styles
   * - handle
     - Pointer to the open DSO file handle

Only one of the three alternate creator entries can be used at a time
and which of those is determined by the style of plugin. The
"creator.v1" element is for factory functions of supported styles
computing forces (i.e.  command, pair, bond, angle, dihedral, or
improper styles) and the function takes as single argument the pointer
to the LAMMPS instance. The factory function is cast to the
``lammpsplugin_factory1`` type before assignment.  The "creator.v2"
element is for factory functions creating an instance of a fix, compute,
or region style and takes three arguments: a pointer to the LAMMPS
instance, an integer with the length of the argument list and a ``char
**`` pointer to the list of arguments. The factory function pointer
needs to be cast to the ``lammpsplugin_factory2`` type before
assignment.

Pair style example
^^^^^^^^^^^^^^^^^^

As an example, a hypothetical pair style plugin "morse2" implemented in
a class ``PairMorse2`` in the files ``pair_morse2.h`` and
``pair_morse2.cpp`` with the factory function and initialization
function would look like this:

.. code-block:: C++

  #include "lammpsplugin.h"
  #include "version.h"
  #include "pair_morse2.h"

  using namespace LAMMPS_NS;

  static Pair *morse2creator(LAMMPS *lmp)
  {
    return new PairMorse2(lmp);
  }

  extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
  {
    lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;
    lammpsplugin_t plugin;

    plugin.version = LAMMPS_VERSION;
    plugin.style   = "pair";
    plugin.name    = "morse2";
    plugin.info    = "Morse2 variant pair style v1.0";
    plugin.author  = "Axel Kohlmeyer (akohlmey@gmail.com)";
    plugin.creator.v1 = (lammpsplugin_factory1 *) &morse2creator;
    plugin.handle  = handle;
    (*register_plugin)(&plugin,lmp);
  }

The factory function in this example is called ``morse2creator()``.  It
receives a pointer to the LAMMPS class as only argument and thus has to
be assigned to the *creator.v1* member of the plugin struct and cast to
the ``lammpsplugin_factory1`` function pointer type.  It returns a
pointer to the allocated class instance derived from the ``Pair`` class.
This function may be declared static to avoid clashes with other
plugins.  The name of the derived class, ``PairMorse2``, however must be
unique inside the entire LAMMPS executable.

Fix style example
^^^^^^^^^^^^^^^^^

If the factory function would be for a fix or compute, which take three
arguments (a pointer to the LAMMPS class, the number of arguments and the
list of argument strings), then the pointer type is ``lammpsplugin_factory2``
and it must be assigned to the *creator.v2* member of the plugin struct.
Below is an example for that:

.. code-block:: C++

  #include "lammpsplugin.h"
  #include "version.h"
  #include "fix_nve2.h"

  using namespace LAMMPS_NS;

  static Fix *nve2creator(LAMMPS *lmp, int argc, char **argv)
  {
    return new FixNVE2(lmp,argc,argv);
  }

  extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
  {
    lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;
    lammpsplugin_t plugin;

    plugin.version = LAMMPS_VERSION;
    plugin.style   = "fix";
    plugin.name    = "nve2";
    plugin.info    = "NVE2 variant fix style v1.0";
    plugin.author  = "Axel Kohlmeyer (akohlmey@gmail.com)";
    plugin.creator.v2 = (lammpsplugin_factory2 *) &nve2creator;
    plugin.handle  = handle;
    (*register_plugin)(&plugin,lmp);
  }

Command style example
^^^^^^^^^^^^^^^^^^^^^
Command styles also use the first variant of factory function as
demonstrated in the following example, which also shows that the
implementation of the plugin class may be within the same source
file as the plugin interface code:

.. code-block:: C++

   #include "lammpsplugin.h"

   #include "comm.h"
   #include "error.h"
   #include "command.h"
   #include "version.h"

   #include <cstring>

   namespace LAMMPS_NS {
     class Hello : public Command {
      public:
       Hello(class LAMMPS *lmp) : Command(lmp) {};
       void command(int, char **);
     };
   }

   using namespace LAMMPS_NS;

   void Hello::command(int argc, char **argv)
   {
      if (argc != 1) error->all(FLERR,"Illegal hello command");
      if (comm->me == 0)
        utils::logmesg(lmp,fmt::format("Hello, {}!\n",argv[0]));
   }

   static void hellocreator(LAMMPS *lmp)
   {
     return new Hello(lmp);
   }

   extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
   {
     lammpsplugin_t plugin;
     lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

     plugin.version = LAMMPS_VERSION;
     plugin.style   = "command";
     plugin.name    = "hello";
     plugin.info    = "Hello world command v1.1";
     plugin.author  = "Axel Kohlmeyer (akohlmey@gmail.com)";
     plugin.creator.v1 = (lammpsplugin_factory1 *) &hellocreator;
     plugin.handle  = handle;
     (*register_plugin)(&plugin,lmp);
   }

Additional Details
^^^^^^^^^^^^^^^^^^

The initialization function **must** be called ``lammpsplugin_init``, it
**must** have C bindings and it takes three void pointers as arguments.
The first is a pointer to the LAMMPS class that calls it and it needs to
be passed to the registration function.  The second argument is a
pointer to the internal handle of the DSO file, this needs to be added
to the plugin info struct, so that the DSO can be closed and unloaded
when all its contained plugins are unloaded.  The third argument is a
function pointer to the registration function and needs to be stored
in a variable of ``lammpsplugin_regfunc`` type and then called with a
pointer to the ``lammpsplugin_t`` struct and the pointer to the LAMMPS
instance as arguments to register a single plugin.  There may be multiple
calls to multiple plugins in the same initialization function.

To register a plugin a struct of the ``lammpsplugin_t`` needs to be filled
with relevant info: current LAMMPS version string, kind of style, name of
style, info string, author string, pointer to factory function, and the
DSO handle.  The registration function is called with a pointer to the address
of this struct and the pointer of the LAMMPS class.  The registration function
will then add the factory function of the plugin style to the respective
style map under the provided name.  It will also make a copy of the struct
in a list of all loaded plugins and update the reference counter for loaded
plugins from this specific DSO file.

The pair style itself (i.e. the PairMorse2 class in this example) can be
written just like any other pair style that is included in LAMMPS.  For
a plugin, the use of the ``PairStyle`` macro in the section encapsulated
by ``#ifdef PAIR_CLASS`` is not needed, since the mapping of the class
name to the style name is done by the plugin registration function with
the information from the ``lammpsplugin_t`` struct.  It may be included
in case the new code is intended to be later included in LAMMPS directly.
