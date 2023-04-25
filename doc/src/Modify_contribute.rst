Submitting new features for inclusion in LAMMPS
===============================================

We encourage LAMMPS users to submit new features they wrote for LAMMPS
to be included into the LAMMPS distribution and thus become easily
accessible to all LAMMPS users.  The LAMMPS source code is managed with
git and public development is hosted on `GitHub
<https://github.com/lammps/lammps>`_.  You can monitor the repository to
be notified of releases, follow the ongoing development, and comment on
topics of interest to you.

Communication with the LAMMPS developers
----------------------------------------

For any larger modifications or programming project, you are encouraged
to contact the LAMMPS developers ahead of time in order to discuss
implementation strategies and coding guidelines.  That will make it
easier to integrate your contribution and results in less work for
everybody involved.  You are also encouraged to search through the list
of `open issues on GitHub <https://github.com/lammps/lammps/issues>`_
and submit a new issue for a planned feature, so you would not duplicate
the work of others (and possibly get scooped by them) or have your work
duplicated by others.

For informal communication with the LAMMPS developers you may ask to
join the `LAMMPS developers on Slack <https://lammps.slack.com>`_.  This
slack work space is by invitation only.  Thus for access, please send an
e-mail to ``slack@lammps.org`` explaining what part of LAMMPS you are
working on.  Only discussions related to LAMMPS development are
tolerated in that work space, so this is **NOT** for people that look
for help with compiling, installing, or using LAMMPS.  Please post a
message to the `LAMMPS forum <https://www.lammps.org/forum.html>`_ for
those purposes.

Packages versus individual files
--------------------------------

The remainder of this chapter describes how to add new "style" files of
various kinds to LAMMPS.  Packages are simply collections of one or more
such new class files which are invoked as a new style within a LAMMPS
input script.  In some cases collections of supporting functions or
classes are also included as separate files in a package, especially when
they can be shared between multiple styles. If designed correctly, these
additions typically do not require any changes to the core code of
LAMMPS; they are simply add-on files that are compiled with the rest of
LAMMPS.  To make those styles work, you may need some trivial changes to
the core code; an example of a trivial change is making a parent-class
method "virtual" when you derive a new child class from it.

If you think your new feature or package requires some non-trivial
changes in core LAMMPS files, you should communicate with the LAMMPS
developers `on Slack <https://lammps.org/slack.html>`_, `on GitHub
<https://github.com/lammps/lammps/issues>`_, or `via email
<https://www.lammps.org/authors.html>`_, since we may have
recommendations about what changes to do where, or may not want to
include certain changes for some reason and thus you would need to look
for alternatives.

Time and effort required
------------------------

How quickly your contribution will be integrated can vary a lot.  It
depends largely on how much effort it will cause the LAMMPS developers
to integrate and test it, how many and what kind of changes to the core
code are required, how quickly you can address them and of how much
interest it is to the larger LAMMPS community.  Please see the section
on :doc:`LAMMPS programming style and requirements <Modify_style>` for
instructions, recommendations, and formal requirements.  A small,
modular, well written contribution may be integrated within hours, but a
complex change that will require a redesign of some core functionality
in LAMMPS for a clean integration can take many months until it is
considered ready for inclusion (though this is rare).


Submission procedure
--------------------

All changes to LAMMPS (including those from LAMMPS developers) are
integrated via pull requests on GitHub and cannot be merged without
passing the automated testing and an approving review by a LAMMPS core
developer.  Thus before submitting your contribution, you should first
make certain, that your added or modified code compiles and works
correctly with the latest development version of LAMMPS and contains all
bug fixes from it.

Once you have prepared everything, see the :doc:`LAMMPS GitHub Tutorial
<Howto_github>` page for instructions on how to submit your changes or
new files through a GitHub pull request yourself.  If you are unable or
unwilling to submit via GitHub yourself, you may also submit patch files
or full files to the LAMMPS developers and ask them to submit a pull
request on GitHub on your behalf.  If this is the case, create a gzipped
tar file of all new or changed files or a corresponding patch file using
'diff -u' or 'diff -c' format and compress it with gzip.  Please only
use gzip compression, as this works well and is available on all platforms.

If the new features/files are broadly useful we may add them as core
files to LAMMPS or as part of a :doc:`package <Packages_list>`.  All
packages are listed and described on the :doc:`Packages details
<Packages_details>` doc page.

Licensing
---------

Note that by providing us files to release, you agree to make them
open-source, i.e. we can release them under the terms of the GPL
(version 2) with the rest of LAMMPS.  And similarly as part of a LGPL
(version 2.1) distribution of LAMMPS that we make available to
developers on request only and with files that are not authorized for
that kind of distribution removed (e.g. interface to FFTW).  See the
:doc:`LAMMPS license <Intro_opensource>` page for details.

External contributions
----------------------

If you prefer to do so, you can also develop and support your add-on
feature **without** having it included in the LAMMPS distribution, for
example as a download from a website of your own.  See the `External
LAMMPS packages and tools <https://www.lammps.org/external.html>`_ page
of the LAMMPS website for examples of groups that do this.  We are happy
to advertise your package and website from that page.  Simply email the
`developers <https://www.lammps.org/authors.html>`_ with info about your
package and we will post it there.  We recommend to name external
packages USER-\<name\> so they can be easily distinguished from bundled
packages that do not have the USER- prefix.

