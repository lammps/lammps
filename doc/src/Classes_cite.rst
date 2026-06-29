Citation management for contributed features
--------------------------------------------

LAMMPS provides a built-in mechanism to remind users to cite publications
that describe the implementation of contributed features. This is
implemented through the :cpp:class:`CiteMe <LAMMPS_NS::CiteMe>` class.

Overview
^^^^^^^^

When users use specific features in LAMMPS, a citation reminder
can be configured in the contributed source code to remind users
which publications they should cite.
Citations are output in three ways:

* To the screen (with configurable verbosity)
* To the log file (with configurable verbosity)
* To an optional BibTeX file for easy integration with bibliography managers

The system automatically deduplicates citations, so each publication is
listed only once even if multiple features reference it or they are used
multiple times.

Adding a citation reminder to contributed code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are contributing a new feature to LAMMPS you can add a reminder
for citing a publication describing the implementation of that feature.
This is a simple two-step process:

**Step 1: Define the citation string**

At the top of your implementation file (typically a ``.cpp`` file),
define a static string containing the citation in BibTeX format. The
string should start with a single line description including a DOI URL,
followed by a complete BibTeX entry. Example:

.. code-block:: c++

   static const char cite_my_feature[] =
     "my_feature command: https://doi.org/10.1234/example.doi\n\n"
     "@Article{AuthorYear,\n"
     " author = {First Author and Second Author},\n"
     " title = {Title of the Paper},\n"
     " journal = {Journal Name},\n"
     " year = 2024,\n"
     " volume = 100,\n"
     " pages = {1-10}\n"
     "}\n\n";

**Step 2: Register the citation**

In your style's constructor, add a call to register the citation. Always
check that the ``citeme`` pointer instance in the LAMMPS is available
before calling, since the pointer will be NULL when citation tracking is
disabled via a command line option:

.. code-block:: c++

   MyStyle::MyStyle(LAMMPS *lmp) : Parent(lmp)
   {
     if (lmp->citeme) lmp->citeme->add(cite_my_feature);
     // ... rest of constructor code
   }

**Output:**

The example from above will produce by default the following output:

.. parsed-literal::

   CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE

   Your simulation uses code contributions which should be cited:

   - my_feature command: https://doi.org/10.1234/example.doi

   The file log.cite lists these citations in BibTeX format.

   CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

The :cpp:class:`CiteMe <LAMMPS_NS::CiteMe>` class uses hash-based
deduplication to ensure each citation is shown only once, even if multiple
features in a simulation reference the same publication.  Only MPI rank 0
performs citation output to avoid flooding the output in parallel runs.

Citations added with :cpp:func:`CiteMe::add() <LAMMPS_NS::CiteMe::add>` are:

1. Written immediately to the BibTeX file (if enabled)
2. Buffered for screen and log file output
3. Flushed to screen and log file at appropriate times (typically at the
   end of the run or when LAMMPS has reached the end of the input file).

For the complete API documentation see the class reference below:

-------------

.. doxygenclass:: LAMMPS_NS::CiteMe
   :project: progguide
   :members:
