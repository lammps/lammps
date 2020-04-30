Parsing rules for input scripts
===============================

Each non-blank line in the input script is treated as a command.
LAMMPS commands are case sensitive.  Command names are lower-case, as
are specified command arguments.  Upper case letters may be used in
file names or user-chosen ID strings.

Here are 6 rules for how each line in the input script is parsed by
LAMMPS:

.. _one:

1. If the last printable character on the line is a "&" character, the
   command is assumed to continue on the next line.  The next line is
   concatenated to the previous line by removing the "&" character and
   line break.  This allows long commands to be continued across two or
   more lines.  See the discussion of triple quotes in :ref:`6 <six>`
   for how to continue a command across multiple line without using "&"
   characters.

.. _two:

2. All characters from the first "#" character onward are treated as
   comment and discarded.  The exception to this rule is described in
   :ref:`6 <six>`. Note that a comment after a trailing "&" character
   will prevent the command from continuing on the next line.  Also note
   that for multi-line commands a single leading "#" will comment out
   the entire command.

   .. code-block:: LAMMPS

      # this is a comment
      timestep  1.0   # this is also a comment

.. _three:

3. The line is searched repeatedly for $ characters, which indicate
   variables that are replaced with a text string.  The exception to
   this rule is described in :ref:`6 <six>`.

   If the $ is followed by text in curly brackets '{}', then the
   variable name is the text inside the curly brackets.  If no curly
   brackets follow the $, then the variable name is the single character
   immediately following the $.  Thus ${myTemp} and $x refer to variables
   named "myTemp" and "x", while "$xx" will be interpreted as a variable
   named "x" followed by an "x" character.

   How the variable is converted to a text string depends on what style
   of variable it is; see the :doc:`variable <variable>` doc page for
   details.  It can be a variable that stores multiple text strings, and
   return one of them.  The returned text string can be multiple "words"
   (space separated) which will then be interpreted as multiple
   arguments in the input command.  The variable can also store a
   numeric formula which will be evaluated and its numeric result
   returned as a string.

   As a special case, if the $ is followed by parenthesis "()", then the
   text inside the parenthesis is treated as an "immediate" variable and
   evaluated as an :doc:`equal-style variable <variable>`.  This is a
   way to use numeric formulas in an input script without having to
   assign them to variable names.  For example, these 3 input script
   lines:

   .. code-block:: LAMMPS

      variable X equal (xlo+xhi)/2+sqrt(v_area)
      region 1 block $X 2 INF INF EDGE EDGE
      variable X delete

   can be replaced by:

   .. code-block:: LAMMPS

      region 1 block $((xlo+xhi)/2+sqrt(v_area)) 2 INF INF EDGE EDGE

   so that you do not have to define (or discard) a temporary variable,
   "X" in this case.

   Additionally, the "immediate" variable expression may be followed by
   a colon, followed by a C-style format string, e.g. ":%f" or ":%.10g".
   The format string must be appropriate for a double-precision
   floating-point value.  The format string is used to output the result
   of the variable expression evaluation.  If a format string is not
   specified a high-precision "%.20g" is used as the default.

   This can be useful for formatting print output to a desired precision:

   .. code-block:: LAMMPS

      print "Final energy per atom: $(pe/atoms:%10.3f) eV/atom"

   Note that neither the curly-bracket or immediate form of variables
   can contain nested $ characters for other variables to substitute
   for.  Thus you may **NOT** do this:

   .. code-block:: LAMMPS

      variable        a equal 2
      variable        b2 equal 4
      print           "B2 = ${b$a}"

   Nor can you specify an expression like "$($x-1.0)" for an immediate
   variable, but you could use $(v_x-1.0), since the latter is valid
   syntax for an :doc:`equal-style variable <variable>`.

   See the :doc:`variable <variable>` command for more details of how
   strings are assigned to variables and evaluated, and how they can
   be used in input script commands.

.. _four:

4. The line is broken into "words" separated by white-space (tabs,
   spaces).  Note that words can thus contain letters, digits,
   underscores, or punctuation characters.

.. _five:

5. The first word is the command name.  All successive words in the line
   are arguments.

.. _six:

6. If you want text with spaces to be treated as a single argument, it
   can be enclosed in either single or double or triple quotes.  A long
   single argument enclosed in single or double quotes can span multiple
   lines if the "&" character is used, as described above.  When the
   lines are concatenated together (and the "&" characters and line
   breaks removed), the text will become a single line.  If you want
   multiple lines of an argument to retain their line breaks, the text
   can be enclosed in triple quotes, in which case "&" characters are
   not needed.  For example:

   .. code-block:: LAMMPS

      print "Volume = $v"
      print 'Volume = $v'
      if "${steps} > 1000" then quit
      variable a string "red green blue &
                      purple orange cyan"
      print """
      System volume = $v
      System temperature = $t
      """

   In each case, the single, double, or triple quotes are removed when
   the single argument they enclose is stored internally.

   See the :doc:`dump modify format <dump_modify>`, :doc:`print
   <print>`, :doc:`if <if>`, and :doc:`python <python>` commands for
   examples.

   A "#" or "$" character that is between quotes will not be treated as
   a comment indicator in :ref:`2 <two>` or substituted for as a
   variable in :ref:`3 <three>`.

.. note::

   If the argument is itself a command that requires a quoted
   argument (e.g. using a :doc:`print <print>` command as part of an
   :doc:`if <if>` or :doc:`run every <run>` command), then single, double, or
   triple quotes can be nested in the usual manner.  See the doc pages
   for those commands for examples.  Only one of level of nesting is
   allowed, but that should be sufficient for most use cases.
