Thermodynamic output options
============================

The ``Thermo`` class computes and prints thermodynamic information to
the screen and log file; see the files ``thermo.cpp`` and ``thermo.h``.

There are four styles defined in ``thermo.cpp``: "one", "multi", "yaml",
and "custom".  The "custom" style allows the user to explicitly list
keywords for individual quantities to print when thermodynamic output is
generated.  The others have a fixed list of keywords.  See the
:doc:`thermo_style <thermo_style>` command for a list of available
quantities.  The formatting of the "custom" style defaults to the "one"
style, but can be adapted using :doc:`thermo_modify line <thermo_modify>`.

The thermo styles (one, multi, etc) are defined by lists of keywords
with associated formats for integer and floating point numbers and
identified but an enumerator constant.  Adding a new style thus mostly
requires defining a new list of keywords and the associated formats and
then inserting the required output processing where the enumerators are
identified.  Search for the word "CUSTOMIZATION" with references to
"thermo style" in the ``thermo.cpp`` file to see the locations where
code will need to be added.  The member function ``Thermo::header()``
prints output at the very beginning of a thermodynamic output block and
can be used to print column headers or other front matter.  The member
function ``Thermo::footer()`` prints output at the end of a
thermodynamic output block.  The formatting of the output is done by
assembling a "line" (which may span multiple lines if the style inserts
newline characters ("\n" as in the "multi" style).

New thermodynamic keywords can also be added to ``thermo.cpp`` to
compute new quantities for output.  Search for the word "CUSTOMIZATION"
with references to "keyword" in ``thermo.cpp`` to see the several
locations where code will need to be added.  Effectively, you need to
define a member function that computes the property, add an if statement
in ``Thermo::parse_fields()`` where the corresponding header string for
the keyword and the function pointer is registered by calling the
``Thermo::addfield()`` method, and add an if statement in
``Thermo::evaluate_keyword()`` which is called from the ``Variable``
class when a thermo keyword is encountered.

.. note::

   The third argument to ``Thermo::addfield()`` is a flag indicating
   whether the function for the keyword computes a floating point
   (FLOAT), regular integer (INT), or big integer (BIGINT) value.  This
   information is used for formatting the thermodynamic output.  Inside
   the function the result must then be stored either in the ``dvalue``,
   ``ivalue`` or ``bivalue`` member variable, respectively.

Since the :doc:`thermo_style custom <thermo_style>` command allows to
use output of quantities calculated by :doc:`fixes <fix>`,
:doc:`computes <compute>`, and :doc:`variables <variable>`, it may often
be simpler to compute what you wish via one of those constructs, rather
than by adding a new keyword to the thermo_style command.
