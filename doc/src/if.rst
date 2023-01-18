.. index:: if

if command
==========

Syntax
""""""

.. code-block:: LAMMPS

   if boolean then t1 t2 ... elif boolean f1 f2 ... elif boolean f1 f2 ... else e1 e2 ...

* boolean = a Boolean expression evaluated as TRUE or FALSE (see below)
* then = required word
* t1,t2,...,tN = one or more LAMMPS commands to execute if condition is met, each enclosed in quotes
* elif = optional word, can appear multiple times
* f1,f2,...,fN = one or more LAMMPS commands to execute if elif condition is met, each enclosed in quotes (optional arguments)
* else = optional argument
* e1,e2,...,eN = one or more LAMMPS commands to execute if no condition is met, each enclosed in quotes (optional arguments)

Examples
""""""""

.. code-block:: LAMMPS

   if "${steps} > 1000" then quit
   if "${myString} == a10" then quit
   if "$x <= $y" then "print 'X is smaller = $x'" else "print 'Y is smaller = $y'"
   if "(${eng} > 0.0) || ($n < 1000)" then &
     "timestep 0.005" &
   elif $n<10000 &
     "timestep 0.01" &
   else &
     "timestep 0.02" &
     "print 'Max step reached'"
   if "${eng} > ${eng_previous}" then "jump file1" else "jump file2"

Description
"""""""""""

This command provides an if-then-else capability within an input
script.  A Boolean expression is evaluated and the result is TRUE or
FALSE.  Note that as in the examples above, the expression can contain
variables, as defined by the :doc:`variable <variable>` command, which
will be evaluated as part of the expression.  Thus a user-defined
formula that reflects the current state of the simulation can be used
to issue one or more new commands.

If the result of the Boolean expression is TRUE, then one or more
commands (t1, t2, ..., tN) are executed.  If it is FALSE, then Boolean
expressions associated with successive elif keywords are evaluated
until one is found to be true, in which case its commands (f1, f2,
..., fN) are executed.  If no Boolean expression is TRUE, then the
commands associated with the else keyword, namely (e1, e2, ..., eN),
are executed.  The elif and else keywords and their associated
commands are optional.  If they are not specified and the initial
Boolean expression is FALSE, then no commands are executed.

The syntax for Boolean expressions is described below.

Each command (t1, f1, e1, etc.) can be any valid LAMMPS input script
command.  If the command is more than one word, it must enclosed in
quotes, so it will be treated as a single argument, as in the examples
above.

.. note::

   If a command itself requires a quoted argument (e.g., a
   :doc:`print <print>` command), then double and single quotes can be used
   and nested in the usual manner, as in the examples above and below.
   The :doc:`Commands parse <Commands_parse>` page has more details on
   using quotes in arguments.  Only one of level of nesting is allowed,
   but that should be sufficient for most use cases.

Note that by using the line continuation character "&", the if command
can be spread across many lines, though it is still a single command:

.. code-block:: LAMMPS

   if "$a < $b" then &
     "print 'Minimum value = $a'" &
     "run 1000" &
   else &
     'print "Minimum value = $b"' &
     "minimize 0.001 0.001 1000 10000"

Note that if one of the commands to execute is :doc:`quit <quit>`, as in
the first example above, then executing the command will cause LAMMPS
to halt.

Note that by jumping to a label in the same input script, the if
command can be used to break out of a loop.  See the :doc:`variable delete <variable>` command for info on how to delete the associated
loop variable, so that it can be re-used later in the input script.

Here is an example of a loop which checks every 1000 steps if the
system temperature has reached a certain value, and if so, breaks out
of the loop to finish the run.  Note that any variable could be
checked, so long as it is current on the timestep when the run
completes.  As explained on the :doc:`variable <variable>` doc page,
this can be insured by including the variable in thermodynamic output.

.. code-block:: LAMMPS

   variable myTemp equal temp
   label loop
   variable a loop 1000
   run 1000
   if "${myTemp} < 300.0" then "jump SELF break"
   next a
   jump SELF loop
   label break
   print "ALL DONE"

Here is an example of a double loop which uses the if and
:doc:`jump <jump>` commands to break out of the inner loop when a
condition is met, then continues iterating through the outer loop.

.. code-block:: LAMMPS

   label       loopa
   variable    a loop 5
     label     loopb
     variable  b loop 5
       print     "A,B = $a,$b"
       run       10000
       if        "$b > 2" then "jump SELF break"
     next      b
     jump      in.script loopb
     label       break
     variable    b delete
   next        a
   jump        SELF loopa

----------

The Boolean expressions for the if and elif keywords have a C-like
syntax.  Note that each expression is a single argument within the if
command.  Thus if you want to include spaces in the expression for
clarity, you must enclose the entire expression in quotes.

An expression is built out of numbers (which start with a digit or
period or minus sign) or strings (which start with a letter and can
contain alphanumeric characters, underscores, or forward slashes):

.. parsed-literal::

   0.2, 100, 1.0e20, -15.4, ...
   InP, myString, a123, ab_23_cd, lj/cut, ...

and Boolean operators:

.. parsed-literal::

   A == B, A != B, A < B, A <= B, A > B, A >= B, A && B, A \|\| B, A \|\^ B, !A

Each A and B is a number or string or a variable reference like ``$a`` or
``${abc}``, or A or B can be another Boolean expression.

Note that all variables used will be substituted for before the
Boolean expression in evaluated.  A variable can produce a number,
like an :doc:`equal-style variable <variable>`, or it can produce a
string, like an :doc:`index-style variable <variable>`.

The Boolean operators ``==`` and ``!=`` can operate on a pair or strings
or numbers.  They cannot compare a number to a string.  All the other
Boolean operations can only operate on numbers.

Expressions are evaluated left to right and have the usual C-style
precedence: the unary logical NOT operator ``!`` has the highest
precedence, the 4 relational operators ``<``, ``<=``, ``>``, and ``>=`` are
next; the two remaining relational operators ``==`` and ``!=`` are next;
then the logical AND operator ``&&``; and finally the logical OR
operator ``||`` and logical XOR (exclusive or) operator ``|^`` have
the lowest precedence.  Parenthesis can be used to group one or more
portions of an expression and/or enforce a different order of
evaluation than what would occur with the default precedence.

When the six relational operators (first six in list above) compare two
numbers, they return either a 1.0 or 0.0 depending on whether the
relationship between A and B is TRUE or FALSE.

When the three logical operators (last three in list above) compare two numbers,
they also return either a 1.0 or 0.0 depending on whether the
relationship between A and B is TRUE or FALSE (or just A).  The
logical AND operator will return 1.0 if both its arguments are
non-zero, else it returns 0.0.  The logical OR operator will return
1.0 if either of its arguments is non-zero, else it returns 0.0.  The
logical XOR operator will return 1.0 if one of its arguments is zero
and the other non-zero, else it returns 0.0.  The logical NOT operator
returns 1.0 if its argument is 0.0, else it returns 0.0.  The 3
logical operators can only be used to operate on numbers, not on
strings.

The overall Boolean expression produces a TRUE result if the numeric
result is non-zero.  If the result is zero, the expression result is
FALSE.

.. note::

   If the Boolean expression is a single numeric value with no Boolean
   operators, it will be FALSE if the value = 0.0, otherwise TRUE.  If
   the Boolean expression is a single string, an error message will be
   issued.

----------

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`variable <variable>`, :doc:`print <print>`

Default
"""""""

none
