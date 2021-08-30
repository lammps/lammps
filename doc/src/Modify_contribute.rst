Submitting new features for inclusion in LAMMPS
===============================================

We encourage users to submit new features or modifications for LAMMPS to
`the core developers <https://www.lammps.org/authors.html>`_ so they can
be added to the LAMMPS distribution. The preferred way to manage and
coordinate this is via the LAMMPS project on `GitHub
<https://github.com/lammps/lammps>`_.  Please see the :doc:`GitHub
Tutorial <Howto_github>` for a demonstration on how to do that.  An
alternative is to contact the LAMMPS developers or the indicated
developer of a package or feature directly and send in your contribution
via e-mail, but that can add a significant delay on getting your
contribution included, depending on how busy the respective developer
is, how complex a task it would be to integrate that code, and how
many - if any - changes are required before the code can be included.

For any larger modifications or programming project, you are encouraged
to contact the LAMMPS developers ahead of time in order to discuss
implementation strategies and coding guidelines. That will make it
easier to integrate your contribution and results in less work for
everybody involved.  You are also encouraged to search through the list
of `open issues on GitHub <https://github.com/lammps/lammps/issues>`_
and submit a new issue for a planned feature, so you would not duplicate
the work of others (and possibly get scooped by them) or have your work
duplicated by others.

For informal communication with the LAMMPS developers you may ask to
join the `LAMMPS developers on Slack <https://lammps.slack.com>`_.  This
slack work space is by invitation only. Thus for access, please send an
e-mail to ``slack@lammps.org`` explaining what part of LAMMPS you are
working on.  Only discussions related to LAMMPS development are
tolerated, so this is **NOT** for people that look for help with
compiling, installing, or using LAMMPS. Please post a message to the
`lammps-users mailing list <https://www.lammps.org/mail.html>`_ or the
`LAMMPS forum <https://www.lammps.org/forum.html>`_ for those purposes.

How quickly your contribution will be integrated depends largely on how
much effort it will cause to integrate and test it, how many and what
kind of changes it requires to the core codebase, and of how much
interest it is to the larger LAMMPS community.  Please see the section
on :doc:`LAMMPS programming style and requirements <Modify_style>` for
instructions, recommendations, and requirements.

Once you have prepared everything, see the :doc:`LAMMPS GitHub Tutorial
<Howto_github>` page for instructions on how to submit your changes or
new files through a GitHub pull request.  If you are unable or unwilling
to submit via GitHub yourself, you may also submit patch files or full
files to the LAMMPS developers and ask them to submit a pull request on
GitHub on your behalf.  **All** changes to LAMMPS (including those from
the LAMMPS core developers) must be submitted as GitHub pull requests
and cannot be merged without passing the automated integration and unit
testing as well as a code review by a LAMMPS core developer that did not
submit it.  Thus before submitting your contribution, you should first
make certain, that your added or modified code works correctly with the
latest patch-level version of LAMMPS and contains all bug fixes from it.
Then create a gzipped tar file of all changed or added files or a
corresponding patch file using 'diff -u' or 'diff -c' and compress it
with gzip.  Please only use gzip compression, as this works well and is
available on all platforms.

If the new features/files are broadly useful we may add them as core
files to LAMMPS or as part of a :doc:`package <Packages_list>`.  All
packages are listed and described on the :doc:`Packages details
<Packages_details>` doc page.

Note that by providing us files to release, you agree to make them
open-source, i.e. we can release them under the terms of the GPL
(version 2) with the rest of LAMMPS.  And similarly as part of a LGPL
(version 2.1) distribution of LAMMPS that we make available to
developers on request only and with files that are not authorized for
that kind of distribution removed (e.g. interface to FFTW).  See the
:doc:`LAMMPS license <Intro_opensource>` page for details.

.. note::

   If you prefer to do so, you can also develop and support your add-on
   feature without having it included in the LAMMPS distribution, for
   example as a download from a website of your own.  See the `Offsite
   LAMMPS packages and tools <https://www.lammps.org/offsite.html>`_
   page of the LAMMPS website for examples of groups that do this.  We
   are happy to advertise your package and website from that page.
   Simply email the `developers <https://www.lammps.org/authors.html>`_
   with info about your package and we will post it there.  We recommend
   to name external packages USER-\<name\> so they can be easily
   distinguished from bundled packages that do not have the USER-
   prefix.

.. _lws: https://www.lammps.org

The previous sections of this page describe how to add new "style"
files of various kinds to LAMMPS.  Packages are simply collections of
one or more new class files which are invoked as a new style within a
LAMMPS input script.  If designed correctly, these additions typically
do not require changes to the main core of LAMMPS; they are simply
add-on files.  If you think your new feature requires non-trivial
changes in core LAMMPS files, you should `communicate with the
developers <https://www.lammps.org/authors.html>`_, since we may or
may not want to include those changes for some reason.  An example of a
trivial change is making a parent-class method "virtual" when you derive
a new child class from it.
