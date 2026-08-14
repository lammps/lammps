Dump styles
===========

Classes that dump per-atom info to files are derived from the Dump
class.  To dump new quantities or in a new format, a new derived dump
class can be added, but it is typically simpler to modify the
DumpCustom class contained in the ``dump_custom.cpp`` file.

``src/dump_atom.cpp`` is a simple example of a derived dump class.

Here is a brief description of methods you define in your new derived
class.  See ``src/dump.h`` for details.

+---------------------+---------------------------------------------------+
| Required            | "pure" methods that *must* be overridden          |
+=====================+===================================================+
| init_style          | initialize style-specific state before a run      |
+---------------------+---------------------------------------------------+
| write_header        | write the header section of a snapshot of atoms   |
+---------------------+---------------------------------------------------+
| pack                | pack a proc's output data into a buffer           |
+---------------------+---------------------------------------------------+
| write_data          | write a proc's data to a file                     |
+---------------------+---------------------------------------------------+

+---------------------+-----------------------------------------------------+
| Optional            | methods that have a default or empty implementation |
+=====================+=====================================================+
| count               | count the number of lines a processor will output   |
+---------------------+-----------------------------------------------------+
| openfile            | open the output file (override for custom I/O)      |
+---------------------+-----------------------------------------------------+
| convert_string      | convert atom data to a string buffer                |
+---------------------+-----------------------------------------------------+
| write_footer        | write a footer after each snapshot                  |
+---------------------+-----------------------------------------------------+
| pack_forward_comm   | pack a buffer for forward ghost communication       |
+---------------------+-----------------------------------------------------+
| unpack_forward_comm | unpack a forward communication buffer               |
+---------------------+-----------------------------------------------------+
| pack_reverse_comm   | pack a buffer for reverse ghost communication       |
+---------------------+-----------------------------------------------------+
| unpack_reverse_comm | unpack a reverse communication buffer               |
+---------------------+-----------------------------------------------------+
| modify_param        | called when dump_modify is executed (optional)      |
+---------------------+-----------------------------------------------------+
| extract             | provide access to internal data (optional)          |
+---------------------+-----------------------------------------------------+

See the :doc:`dump <dump>` command and its *custom* style for a list of
keywords for atom information that can already be dumped by
DumpCustom.  It includes options to dump per-atom info from Compute
classes, so adding a new derived Compute class is one way to calculate
new quantities to dump.

Note that new keywords for atom properties are not typically
added to the :doc:`dump custom <dump>` command.  Instead they are added
to the :doc:`compute property/atom <compute_property_atom>` command.
