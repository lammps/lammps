.. index:: print

print command
=============

Syntax
""""""


.. parsed-literal::

   print string keyword value

* string = text string to print, which may contain variables
* zero or more keyword/value pairs may be appended
* keyword = *file* or *append* or *screen* or *universe*

  .. parsed-literal::

       *file* value = filename
       *append* value = filename
       *screen* value = *yes* or *no*
       *universe* value = *yes* or *no*



Examples
""""""""


.. parsed-literal::

   print "Done with equilibration" file info.dat
   print Vol=$v append info.dat screen no
   print "The system volume is now $v"
   print 'The system volume is now $v'
   print "NEB calculation 1 complete" screen no universe yes
   print """
   System volume = $v
   System temperature = $t
   """

Description
"""""""""""

Print a text string to the screen and logfile.  The text string must
be a single argument, so if it is one line but more than one word, it
should be enclosed in single or double quotes.  To generate multiple
lines of output, the string can be enclosed in triple quotes, as in
the last example above.  If the text string contains variables, they
will be evaluated and their current values printed.

If the *file* or *append* keyword is used, a filename is specified to
which the output will be written.  If *file* is used, then the
filename is overwritten if it already exists.  If *append* is used,
then the filename is appended to if it already exists, or created if
it does not exist.

If the *screen* keyword is used, output to the screen and logfile can
be turned on or off as desired.

If the *universe* keyword is used, output to the global screen and
logfile can be turned on or off as desired. In multi-partition
calculations, the *screen* option and the corresponding output only
apply to the screen and logfile of the individual partition.

If you want the print command to be executed multiple times (with
changing variable values), there are 3 options.  First, consider using
the :doc:`fix print <fix_print>` command, which will print a string
periodically during a simulation.  Second, the print command can be
used as an argument to the *every* option of the :doc:`run <run>`
command.  Third, the print command could appear in a section of the
input script that is looped over (see the :doc:`jump <jump>` and
:doc:`next <next>` commands).

See the :doc:`variable <variable>` command for a description of *equal*
style variables which are typically the most useful ones to use with
the print command.  Equal-style variables can calculate formulas
involving mathematical operations, atom properties, group properties,
thermodynamic properties, global values calculated by a
:doc:`compute <compute>` or :doc:`fix <fix>`, or references to other
:doc:`variables <variable>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix print <fix_print>`, :doc:`variable <variable>`

Default
"""""""

The option defaults are no file output, screen = yes, and universe = no.
