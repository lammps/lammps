Submitting new features for inclusion in LAMMPS
===============================================

We encourage LAMMPS users to submit new features they wrote for LAMMPS
to be included in the LAMMPS distribution and thus become easily
accessible to all LAMMPS users.  The LAMMPS source code is managed
with git and public development is hosted on `GitHub
<https://github.com/lammps/lammps>`_.  You can monitor the repository
to be notified of releases, follow the ongoing development, and
comment on topics of interest to you.

This section contains general information regarding the preparation
and submission of new features to LAMMPS. If you are new to
development in LAMMPS, we recommend you read one of the tutorials on
developing a new :doc:`pair style <Developer_write_pair>` or :doc:`fix
style <Developer_write_fix>` which provide a friendly introduction to
what LAMMPS development entails and common vocabulary used on this
section.


Communication with the LAMMPS developers
----------------------------------------

For any larger modifications or programming project, you are
encouraged to contact the LAMMPS developers ahead of time to discuss
implementation strategies. That will make it easier to integrate your
contribution and typically results in less work for everyone involved.
You are also encouraged to search through the list of `open issues on
GitHub <https://github.com/lammps/lammps/issues>`_ and submit a new
issue for a planned feature, to avoid duplicating work (and possibly
being scooped).

For informal communication with the LAMMPS developers, you may ask to
join the `LAMMPS developers on Slack <https://lammps.slack.com>`_.
This slack work space is by invitation only.  For access, please send
an e-mail to ``slack@lammps.org`` explaining what part of LAMMPS you
are working on.  Only discussions related to LAMMPS development are
tolerated in that work space, so this is **NOT** for people looking
for help with compiling, installing, or using LAMMPS.  Please post a
message to the `LAMMPS forum <https://www.lammps.org/forum.html>`_ for
those purposes.


Time and effort required
------------------------

How quickly your contribution will be integrated can vary widely.  It
depends largely on how much effort is required by the LAMMPS
developers to integrate and test it, if any and what kind of changes
to the core code are required, how quickly you can address them, and
how much interest the contribution is to the larger LAMMPS
community. This process can be streamlined by following the
:doc:`requirements <Modify_requirements>` and :doc:`style
guidelines<Modify_style>`.  A small, modular, well written
contribution may be integrated within hours, but a complex change that
requires a re-design of a core functionality in LAMMPS can take months
before inclusion (though this is rare).


Submission procedure
--------------------

All changes to LAMMPS (including those from LAMMPS developers) are
integrated via pull requests on GitHub and cannot be merged without
passing the automated testing and an approving review by a LAMMPS core
developer.  Before submitting your contribution, you should therefore
first ensure that your added or modified code compiles and works
correctly with the latest development version of LAMMPS and contains
all bug fixes from it.

Once you have prepared everything, see the :doc:`LAMMPS GitHub
Tutorial <Howto_github>` page for instructions on how to submit your
changes or new files through a GitHub pull request.  If you are unable
or unwilling to submit via GitHub yourself, you may also send patch
files or full files to the `LAMMPS developers
<https://www.lammps.org/authors.html>`_ and ask them to submit a pull
request on GitHub on your behalf.  If this is the case, create a
gzipped tar file of all new or changed files or a corresponding patch
file using 'diff -u' or 'diff -c' format and compress it with gzip.
Please only use gzip compression, as this works well and is available
on all platforms.  This mode of submission may delay the integration
as it depends more on the LAMMPS developers.


External contributions
----------------------

If you prefer to do so, you can also develop and support your add-on
feature **without** having it included in the LAMMPS distribution, for
example as a download from a website of your own.  See the `External
LAMMPS packages and tools <https://www.lammps.org/external.html>`_
page of the LAMMPS website for examples of groups that do this.  We
are happy to advertise your package and website from that page.
Simply email the `developers <https://www.lammps.org/authors.html>`_
with info about your package, and we will post it there.  We recommend
naming external packages USER-\<name\> so they can be easily
distinguished from packages in the LAMMPS distribution which do not
have the USER- prefix.


Location of files: individual files and packages
------------------------------------------------

We rarely accept new styles in the core src folder.  Thus, please
review the list of :doc:`available Packages <Packages_details>` to see
if your contribution should be added to one of them.  It should fit
into the general purpose of that package.  If it does not fit well, it
may be added to one of the EXTRA- packages or the MISC package.

However, if your project includes many related features that are not
covered by one of the existing packages or is dependent on a library
(bundled or external), it is best to create a new package with its own
directory (with a name like FOO).  In addition to your new files, the
directory should contain a README text file containing your name and
contact information and a brief description of what your new package
does.


Changes to core LAMMPS files
--------------------------------

If designed correctly, most additions do not require any changes to
the core code of LAMMPS; they are simply add-on files that are
compiled with the rest of LAMMPS.  To make those styles work, you may
need some trivial changes to the core code.  An example of a trivial
change is making a parent-class method "virtual" when you derive a new
child class from it.  If your features involve more substantive
changes to the core LAMMPS files, it is particularly encouraged that
you communicate with the LAMMPS developers early in development.
