Writing or updating documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every contribution that adds a new style or command, or that changes the
behavior or the arguments of an existing one, **must** come with matching
documentation.  This is a :ref:`strict requirement <ReqDocumentation>`
for inclusion into the LAMMPS distribution.  The documentation sources
are written in `reStructuredText <rst_>`_ (files with a ``.rst`` suffix
in the ``doc/src/`` folder) and translated to HTML and PDF with `Sphinx
<https://www.sphinx-doc.org>`_.  All text must be written in American
English and the ``.rst`` files may only contain 7-bit ASCII characters,
so they can be cleanly converted to PDF via PDFLaTeX.  Special characters
and formulas can be included as embedded math typeset in a LaTeX subset.

.. _rst: https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html

This page complements the previous sections on :doc:`writing new styles
<Developer_write>`: where those describe the *C++ code* of a style, this
one describes the *documentation* that has to accompany it.  We use the
:doc:`pair_style born/gauss <pair_born_gauss>` page (file
``doc/src/pair_born_gauss.rst``) and the :doc:`fix nve <fix_nve>` page
(file ``doc/src/fix_nve.rst``) as running examples, since both follow the
established conventions closely and are short enough to quote in full.

A documentation page is not complete by itself: a new style also has to
be registered in the corresponding command tables and overview lists, and
the build and consistency checks described at the :ref:`bottom of this
page <doc_checks>` must pass without errors or warnings.

Anatomy of a style documentation page
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical per-style documentation page is built from a fixed sequence of
elements.  The order and the exact wording of the section headers matter,
because several of them are verified by automated checks and because
common include files (like ``accel_styles.rst``) are shared across all
pages.  The most reliable way to create a new page is to copy the page of
a similar, existing style and adapt it.

Index entries
"""""""""""""

The very first lines of a documentation page are the Sphinx index
directives.  There must be one ``.. index::`` entry for **every** style
that is documented on the page, *including every accelerated variant*
(``/gpu``, ``/intel``, ``/kk``, ``/omp``, ``/opt``).

.. code-block:: rst

   .. index:: pair_style born/gauss
   .. index:: pair_style born/gauss/kk
   .. index:: pair_style born/gauss/omp

A single documentation file may describe several closely related styles
(for example ``doc/src/pair_born.rst`` documents ``born``,
``born/coul/long``, ``born/coul/msm``, and ``born/coul/wolf``).  In that
case each of those styles, and each of their accelerated variants, needs
its own ``.. index::`` line.  The completeness of these index entries is
verified by ``make style_check`` (see :ref:`below <doc_checks>`), which
parses the ``...Style(name,class)`` macros from the C++ headers and
reports any style or accelerated variant whose index entry is missing.

Title and command heading
"""""""""""""""""""""""""

The page title is the name of the command followed by the word
"command", underlined with a row of ``=`` characters that must be at
least as long as the title text.

.. code-block:: rst

   pair_style born/gauss command
   =============================

When a page documents multiple styles, it is customary to list the
primary style in the title (or to use a more generic title), and to make
sure each documented style has an index entry and is reachable through
the command tables.

In reStructuredText a heading is indicated by underlining it with
characters and the the level of a heading is set by kind of character
used.  The underline must be **at least as long** as the heading text; a
shorter underline will trigger an error message, and the convention is
to match the length exactly.  LAMMPS uses a fixed set of characters for
the heading levels so that pages nest consistently.  This matters in
particular for the combined PDF/LaTeX manual, where all pages are
concatenated into a single document and therefore have to share one
global hierarchy.  From the outermost to the innermost level the
characters are:

* ``======`` -- the title of an individual page
* ``------`` -- a section
* ``^^^^^^`` -- a subsection
* ``""""""`` -- a subsubsection

A per-style reference page needs only two of these levels: ``======``
for the page title (the line ending in "command") and ``""""""`` for
each of the standard section headers (Syntax, Examples, Description,
...).  Pages with a deeper structure -- such as Howto or Programmer
Guide pages, and this developer page -- additionally use ``------`` and
``^^^^^^`` for the intermediate levels.

The page title and the section headers show up as collapsible,
expandable entries in the navigation sidebar of the HTML manual (the
theme is configured with a navigation depth of three levels) and in the
table of contents of the PDF manual.

Accelerator Variants line
"""""""""""""""""""""""""

If, and only if, one or more accelerated variants of the style exist,
an "Accelerator Variants:" line follows the title.  It lists the full
names of the accelerated styles in italics.

.. code-block:: rst

   Accelerator Variants: *born/gauss/kk*, *born/gauss/omp*

This line and the ``.. index::`` entries for the accelerated variants
must be kept consistent with each other and with the actual
accelerated styles present in the source tree.  Further down the page,
the shared file ``accel_styles.rst`` is included (see :ref:`below
<doc_accel_include>`); ``make style_check`` enforces that a page which
lists accelerated index entries also includes ``accel_styles.rst`` and,
conversely, that a page without accelerated styles does not.

Syntax, Examples, and Description
"""""""""""""""""""""""""""""""""

The body of the page begins with a "Syntax" section.  It shows the
command in a ``LAMMPS`` code block and then explains each argument with
a bullet list.

.. code-block:: rst

   Syntax
   """"""

   .. code-block:: LAMMPS

      pair_style born/gauss cutoff

   * born/gauss = name of the pair style
   * cutoff = global cutoff (distance units)

The "Examples" section gives one or more short, realistic command
sequences, again in a ``LAMMPS`` code block.

.. code-block:: rst

   Examples
   """"""""

   .. code-block:: LAMMPS

      pair_style born/gauss 10.0
      pair_coeff 1 1 8.2464e13 12.48 0.042644277 0.44 3.56

The "Description" section is the main explanatory text.  It describes
the model or behavior, typically including the functional form as a
``.. math::`` block and a list of the coefficients with their units.

.. code-block:: rst

   Description
   """""""""""

   .. versionadded:: 4Jul2026

   Pair style *born/gauss* computes pairwise interactions ...

   .. math::

      E = A_0 \exp \left( -\alpha r \right) - A_1 \exp\left[ -\beta \left(r - r_0 \right)^2 \right]

Version markers
"""""""""""""""

New publicly visible commands or newly added keywords must be marked
with ``.. versionadded:: TBD``; behavior that changed for an existing
command uses ``.. versionchanged:: TBD``.  The literal text ``TBD`` is
replaced automatically with the upcoming version date at release time.
These directives are standalone markers: any explanatory text follows
as an ordinary, *unindented* paragraph below the directive, **not** as
the indented body of the directive.  Version markers do not apply to
internal (uppercase) style names.

Accelerator styles include
""""""""""""""""""""""""""

.. _doc_accel_include:

For any page that has accelerated variants, a standard block describing
how to use accelerated styles is pulled in from the shared file
``accel_styles.rst``, set off by horizontal rules:

.. code-block:: rst

   ----------

   .. include:: accel_styles.rst

   ----------

Using the shared include keeps the description of the accelerator
suffixes uniform across all pages, so this text should never be written
out by hand.

Category-specific information sections
""""""""""""""""""""""""""""""""""""""

After the description, most styles have one or more information sections
whose header and content are largely standardized per style category.
Copy the wording from a closely related style and adjust it for the
behavior of your style.  Common examples are:

* For **pair styles**: a section titled "Mixing, shift, table, tail
  correction, restart, rRESPA info" describing whether the style
  supports :doc:`pair_modify <pair_modify>` mixing, the shift, table,
  and tail options, whether it writes to :doc:`binary restart files
  <restart>`, and which :doc:`run_style respa <run_style>` keywords it
  supports.
* For **fix styles**: a section titled "Restart, fix_modify, output, run
  start/stop, minimize info" describing what (if anything) the fix
  writes to restart files, which :doc:`fix_modify <fix_modify>` options
  it honors, what global or per-atom quantities it makes available to
  :doc:`output commands <Howto_output>`, and whether it is invoked
  during :doc:`energy minimization <minimize>`.
* For **compute styles**: an "Output info" section describing the shape
  and meaning of the computed scalar, vector, or array and its units.

These sections should be present even when the answer is "none" or "not
supported"; the :doc:`fix nve <fix_nve>` page is a good template for the
minimal case.

Restrictions, related commands, and defaults
""""""""""""""""""""""""""""""""""""""""""""

Three short sections close out the page.

The "Restrictions" section documents any conditions that must be met for
the style to be usable.  The most common one is a **package
requirement**: the style is only available when LAMMPS was built with a
particular package, with a pointer to the :doc:`Build package
<Build_package>` page.

.. code-block:: rst

   Restrictions
   """"""""""""

   This pair style is only enabled if LAMMPS was built with the EXTRA-PAIR
   package.  See the :doc:`Build package <Build_package>` page for more
   info.

Restrictions are *not* limited to package membership, however.  Many
styles impose additional conditions; copy the exact wording from the
page of a similar style.  The most common categories are:

* **Package requirement:** only enabled when built with a given package
  (almost all styles outside the core); always link to :doc:`Build
  package <Build_package>`.
* **Newton setting:** the style requires the :doc:`newton <newton>`
  setting to be "on" (common for many-body and communication-heavy pair
  styles) or, more rarely, "off".
* **Command ordering:** the command must appear before or after another
  command -- for example before the simulation box is defined by
  :doc:`create_box <create_box>` or :doc:`read_data <read_data>`, after
  a prerequisite fix, or before any :doc:`pair_coeff <pair_coeff>`
  commands.
* **Atom style requirement:** the style needs an :doc:`atom style
  <atom_style>` that stores particular per-atom data (e.g. charge,
  finite size/radius, torque and angular momentum, density or
  temperature, or support for bonds and angles).
* **Dependency on another command:** the style only works in conjunction
  with a specific fix, pair style, bond style, or
  :doc:`kspace_style <kspace_style>`.
* **Uniqueness or mutual exclusivity:** only a single instance of the
  style may be defined, or certain keyword combinations cannot be used
  together.
* **Geometry and dimensionality:** the style only works for 2d or only
  for 3d systems, or only for orthogonal (not triclinic) simulation
  boxes.
* **Integrator and parallel constraints:** the style is incompatible
  with :doc:`run_style respa <run_style>`, is not invoked during
  :doc:`energy minimization <minimize>`, or has limitations under MPI
  domain decomposition.

If there are no restrictions at all, the section contains only the
single word ``none``.

The "Related commands" section is a comma-separated list of links
(``:doc:`` roles) to closely related commands.  The "Default" section
describes the default settings of any optional keywords, or the single
word ``none`` if there are none.

Typical section order
"""""""""""""""""""""

Putting it together, the canonical order of a per-style page is:

#. ``.. index::`` entries (one per documented style and accelerated variant)
#. title with the ``command`` suffix
#. "Accelerator Variants:" line (only if accelerated styles exist)
#. Syntax
#. Examples
#. Description (with a ``versionadded``/``versionchanged`` marker if applicable)
#. ``.. include:: accel_styles.rst`` block (only if accelerated styles exist)
#. category-specific info section(s) (Mixing... / Restart... / Output info / ...)
#. Restrictions
#. Related commands
#. Default
#. literature citations (if any)

Literature citations
""""""""""""""""""""

If the page cites publications, the citations are placed at the very
bottom of the file behind a horizontal rule, each with a unique anchor
label and a bold key.

.. code-block:: rst

   ----------

   .. _Bomont:

   **(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)

The in-text reference then uses ``:ref:`(Bomont) <Bomont>```.  Citation
anchor labels must be **unique across all** ``.rst`` files, so a
distinctive label (often including a number or the style name) is
needed; ``make anchor_check`` reports duplicate labels.  See the bottom
of ``doc/fix_nh.rst`` for more elaborate examples of citation
formatting.

Integrating a new style into the command lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating the per-style page is not enough: a new style must also be
added to the compact command table and to the alphabetical overview
list for its category, or the documentation will build with warnings and
``make style_check`` will fail.  For a pair style, two files need an
entry:

* ``Commands_pair.rst`` contains a multi-column table.  Add the style in
  alphabetical order; the accelerator suffix letters are appended in
  parentheses (``g`` = GPU, ``i`` = INTEL, ``k`` = KOKKOS, ``o`` =
  OPENMP, ``t`` = OPT) and the link target is the documentation file:

  .. code-block:: rst

     * :doc:`born/gauss (ko) <pair_born_gauss>`

* ``pair_style.rst`` contains the alphabetical overview list with a
  short, one-line description after the link:

  .. code-block:: rst

     * :doc:`born/gauss <pair_born_gauss>` - Born-Mayer / Gaussian potential

The other style categories follow the same pattern with different
files.  ``make style_check`` (via ``utils/check-styles.py``) checks
exactly these files for completeness:

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * - Category
     - Compact table file
     - Overview list file
   * - pair
     - ``Commands_pair.rst``
     - ``pair_style.rst``
   * - bond
     - ``Commands_bond.rst``
     - ``bond_style.rst``
   * - angle
     - ``Commands_bond.rst``
     - ``angle_style.rst``
   * - dihedral
     - ``Commands_bond.rst``
     - ``dihedral_style.rst``
   * - improper
     - ``Commands_bond.rst``
     - ``improper_style.rst``
   * - compute
     - ``Commands_compute.rst``
     - ``compute.rst``
   * - fix
     - ``Commands_fix.rst``
     - ``fix.rst``
   * - kspace
     - ``Commands_kspace.rst``
     - (none)
   * - dump
     - ``Commands_dump.rst``
     - (none)
   * - command
     - ``Commands_all.rst``
     - (none)

The suffix letters in the compact table must be listed in the canonical
order ``g i k o t``, with no separators (for example ``ko``, never
``ok`` or ``o,k``).

Documenting an accelerated variant of an existing style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a new accelerated variant (e.g. a ``/kk`` or ``/omp`` version) is
added for a style that already exists and is already documented, several
small edits are required rather than a new page:

#. Add a ``.. index::`` line for the new accelerated style near the top
   of the existing page.
#. Add or extend the "Accelerator Variants:" line below the title so it
   lists the new variant.
#. Make sure the page includes the shared ``accel_styles.rst`` block; if
   the style previously had no accelerated variants, add that include
   (with surrounding horizontal rules).
#. Add the corresponding suffix **letter** to the style's entry in the
   compact command table (``Commands_pair.rst``, ``Commands_fix.rst``,
   ...), keeping the canonical ``g i k o t`` letter order.

The overview list file (``pair_style.rst``, ``fix.rst``, ...) does not
encode the suffix letters and therefore does not need to change for a
new accelerated variant.

Adding and referencing example inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For most new features it is preferred (see :ref:`the examples
requirement <ReqExamples>`) to add a small, fast example input under the
``examples`` or ``examples/PACKAGES`` directory.  The documentation page
should point the reader to it.  Two conventions are in common use; both
are acceptable:

* mention the directory inside the "Description" or "Examples" section,
  for instance a line like ``Example input scripts available:
  examples/PACKAGES/foo``;
* reference specific input files by their path in the relevant
  discussion.

Other documentation that may need updating
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on the nature of the contribution, additional documentation
may have to be created or updated:

* **New package:** add a description to ``Packages_details.rst`` and an
  entry to the package list, and, if the package needs special build
  steps, add the corresponding build instructions.
* **Complex features:** if a feature needs more background or a
  multi-command workflow to be used correctly, contribute a
  :doc:`Howto document <Howto>` and add it to the toctree in
  ``Howto.rst``.
* **Library or API changes:** additions or changes to the C++ library
  interface or the Fortran module require doxygen comments in the source
  and matching updates to the "Programmer Guide" section of the manual.
* **Error explanations:** errors with multiple possible causes can be
  given a detailed explanation on the ``Error_details`` page with an
  error code anchor, referenced from the code via
  :cpp:func:`utils::errorurl() <LAMMPS_NS::utils::errorurl>` (see the
  :ref:`error message requirements <ReqErrorMessages>`).
* **Figures and PDFs:** inline figures live in ``doc/JPG`` and
  supplementary PDF files in ``doc/PDF``; see existing pages for how to
  embed them.

Building and checking the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _doc_checks:

Before submitting, build and proofread the documentation from within the
``doc`` directory:

.. code-block:: bash

   make html          # build the HTML manual; there should be no warnings
   make spelling      # run the spell checker
   make anchor_check  # check for duplicate anchors and labels
   make style_check   # check completeness of style tables and index entries
   make check         # run all consistency checks at once

The ``make check`` target bundles ``anchor_check``, ``style_check``,
``package_check``, ``char_check``, and ``role_check``.  The ``char_check``
target rejects any non-ASCII character in the ``.rst`` files, and the
``role_check`` target catches malformed roles and directives (for
example a ``:doc:`` or ``:ref:`` role not immediately followed by a
backtick, a ``text <target>`` hyperlink missing its role, or a directive
written with a single colon instead of ``::``).  Carefully read and
proofread the generated HTML page before opening the pull request.

If the spell checker reports a false positive (such as an author surname
or an established acronym), add it to
``doc/utils/sphinx-config/false_positives.txt``.
