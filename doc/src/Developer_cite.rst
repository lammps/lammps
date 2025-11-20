Citation management for contributed features
-------------------------------------------

LAMMPS provides a built-in mechanism to remind users to cite relevant
publications when they use specific contributed features. This is
implemented through the :cpp:class:`CiteMe <LAMMPS_NS::CiteMe>` class,
which manages citation tracking and output formatting.

Overview
^^^^^^^^

When users enable and use contributed packages or special features in
LAMMPS, the citation system automatically tracks which publications should
be cited and displays appropriate reminders during the simulation run.
Citations are output in three ways:

* To the screen (with configurable verbosity)
* To the log file (with configurable verbosity)
* To an optional BibTeX file for easy integration with bibliography managers

The system automatically deduplicates citations, so each publication is
listed only once even if multiple features reference it.

Adding citation reminders to contributed code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are contributing a new feature to LAMMPS that implements a published
method or algorithm, you should add a citation reminder to ensure users
properly acknowledge the scientific work. This is a simple two-step process:

**Step 1: Define the citation string**

At the top of your implementation file (typically a ``.cpp`` file), define a
static string containing the citation in BibTeX format. The string should
start with a brief description including a DOI (Digital Object Identifier),
followed by a complete BibTeX entry:

.. code-block:: c++

   static const char cite_my_feature[] =
     "my_feature command: doi:10.1234/example.doi\n\n"
     "@Article{AuthorYear,\n"
     " author = {First Author and Second Author},\n"
     " title = {Title of the Paper},\n"
     " journal = {Journal Name},\n"
     " year = 2024,\n"
     " volume = 100,\n"
     " pages = {1-10}\n"
     "}\n\n";

.. note::

   The first line should be concise and include the DOI. This line is shown
   in TERSE output mode, while the full BibTeX entry is shown in VERBOSE mode.

**Step 2: Register the citation**

In your style's constructor, add a call to register the citation. Always
check that the ``citeme`` pointer is not NULL before calling:

.. code-block:: c++

   MyStyle::MyStyle(LAMMPS *lmp) : Parent(lmp)
   {
     if (lmp->citeme) lmp->citeme->add(cite_my_feature);
     // ... rest of constructor code
   }

The check ``if (lmp->citeme)`` is important because citation tracking can
be disabled via command-line options.

Example: Adding a citation to a pair style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a complete example from the ASPHERE package showing how the Gay-Berne
pair style adds its citation:

.. code-block:: c++

   // At the top of pair_gayberne.cpp
   static const char cite_pair_gayberne[] =
     "pair gayberne command: doi:10.1063/1.3058435\n\n"
     "@Article{Brown09,\n"
     " author =  {W. M. Brown and M. K. Petersen and S. J. Plimpton and G. S. Grest},\n"
     " title =   {Liquid Crystal Nanodroplets in Solution},\n"
     " journal = {J.~Chem.\\ Phys.},\n"
     " year =    2009,\n"
     " volume =  130,\n"
     " number =  4,\n"
     " pages =   {044901}\n"
     "}\n\n";

   // In the constructor
   PairGayBerne::PairGayBerne(LAMMPS *lmp) : Pair(lmp)
   {
     if (lmp->citeme) lmp->citeme->add(cite_pair_gayberne);
     // ... rest of constructor
   }

When a user runs a simulation using the Gay-Berne pair style, they will see
output like:

.. parsed-literal::

   CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE

   Your simulation uses code contributions which should be cited:

   - pair gayberne command: doi:10.1063/1.3058435

   The file log.cite lists these citations in BibTeX format.

   CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE

Citation output modes
^^^^^^^^^^^^^^^^^^^^^

LAMMPS supports different verbosity levels for citation output:

* **VERBOSE mode**: Displays the full BibTeX entry
* **TERSE mode**: Displays only the first line (typically the DOI and brief description)

The output mode can be controlled separately for screen and log file output
via command-line options when starting LAMMPS. For example:

.. code-block:: bash

   # Terse on screen, verbose in log file (default)
   lmp -cite screen/terse -cite logfile/verbose -cite file citations.bib

   # Verbose everywhere
   lmp -cite screen/verbose -cite logfile/verbose

   # Minimal output (terse everywhere)
   lmp -cite terse

The citation file (if specified with ``-cite file <filename>``) always
contains the full BibTeX entries regardless of the verbosity settings.

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

The :cpp:class:`CiteMe <LAMMPS_NS::CiteMe>` class uses hash-based
deduplication to ensure each citation is shown only once, even if multiple
features in a simulation reference the same publication. Only MPI rank 0
performs citation output to avoid flooding the output in parallel runs.

Citations added with :cpp:func:`CiteMe::add() <LAMMPS_NS::CiteMe::add>` are:

1. Written immediately to the BibTeX file (if enabled)
2. Buffered for screen and log file output
3. Flushed to screen and log file at appropriate times (typically at the
   end of the run or when :cpp:func:`CiteMe::flush()
   <LAMMPS_NS::CiteMe::flush>` is called)

For complete API documentation, see the :cpp:class:`CiteMe
<LAMMPS_NS::CiteMe>` class reference.

Best practices
^^^^^^^^^^^^^^

When adding citation reminders to your contributed code:

* **Always include a DOI** in the first line of the citation string for easy
  lookup
* **Use proper BibTeX formatting** to ensure citations can be imported into
  bibliography managers
* **Place the citation string** at file scope (not inside a namespace) as a
  static const char array
* **Add citations in constructors** where the feature is first activated, not
  in every method call
* **Test your citation** by running a simple example and checking that the
  output is formatted correctly

Remember that proper citation is essential for academic software. It helps:

* Give credit to the developers and researchers who created the methods
* Allow users to find detailed descriptions of the algorithms
* Track the impact and usage of scientific software contributions
* Encourage continued development and sharing of computational methods
