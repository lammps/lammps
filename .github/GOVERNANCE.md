---
title: "Governance"
description: "How the LAMMPS project is managed"
---

All information in this file is pursuant to the LAMMPS Project
[Technical Charter](https://www.lammps.org/about/charter).  In case of
any conflict or difference in language, the LAMMPS Project Technical
Charter takes priority over this document.

Adopted August 18, 2026.

---

## 1. Technical Steering Committee (TSC)

<br>

### 1.1 Role

The role of the Technical Steering Committee (TSC) is to provide
technical direction to the project.  The TSC will vote on any matters on
which the community is unable to reach consensus.

### 1.2 Membership

Members can be added to the TSC by a majority vote of the TSC.  Members
may be removed from the TSC by a 2/3rd majority vote of the TSC.  If a
TSC member has been inactive for over 6 months, the TSC must hold a vote
on whether to remove that member from the TSC.

### 1.3 Current TSC Members

1. Steve Plimpton ([@sjplimp](https://github.com/sjplimp)), unaffiliated
2. Axel Kohlmeyer ([@akohlmey](https://github.com/akohlmey)), Temple University
3. Aidan Thompson ([@athomps](https://github.com/athomps)), Sandia National Laboratories
4. Richard Berger ([@rbberger](https://github.com/rbberger)), Los Alamos National Laboratory
5. Stan Moore ([@stanmoore1](https://github.com/stanmoore1)), Sandia National Laboratories
6. Joel Clemmer, ([@jtclemm](https://github.com/jtclemm)), Sandia National Laboratories
7. Jake Gissinger ([@jrgissing](https://github.com/jrgissing)), Stevens Institute of Technology
8. Larry Fried (), Lawrence Livermore National Laboratory
9. Chris Knight ([@cjknight](https://github.com/cjknight)), Argonne National Laboratory
10. Evan Weinberg ([@weinbe2](https://github.com/weinbe2)), NVIDIA Corporation
11. Ellad Tadmor ([@tadmor](https://github.com/tadmor)), University of Minnesota
12. Tim Mattox ([@timattox](https://github.com/timattox)), Hewlett Packard Enterprise
13. Ian Bogle ([@IanBogle](https://github.com/IanBogle)), AMD

### 1.4 TSC Chair

The TSC will elect a chair.  The TSC chair runs TSC meetings and may
make interim decisions on urgent matters on behalf of the TSC, which may
be reviewed by the TSC at its next meeting.

Current chair: Aidan Thompson ([@athomps](https://github.com/athomps)), SNL

### 1.5 Meetings and Notes

The TSC meetings are held quarterly.  These meetings are open to the
public, and are held virtually.

More information about how to attend is on the LAMMPS homepage at
[Public Meetings](https://www.lammps.org/community/meetings/).

Meeting notes will be taken and can be found in the
[LAMMPS TSC Meeting Notes](https://github.com/lammps/tsc-meeting-notes)
GitHub repository.

---

## 2. Other Public LAMMPS fora

<br>

### 2.1 MatSci.org Discourse Forum

LAMMPS has joined the [MatSci.org Discouse](https://matsci.org/) community
website together with many other related software packages.
There is a [LAMMPS forum](https://matsci.org/lammps)
which is the main place for communication among LAMMPS users and developers.
This also includes a subcategories for beginners, developers, installation,
and an archive of the lammps-users mailing list that was active from 2005
until 2022.

### 2.2 Monthly Developer Meetings

Public LAMMPS developer meetings are held virtually on the second
Tuesday of every month that has not a TSC meeting, at 12 noon Mountain
Time (i.e. local time in Albuquerque, NM).

More information about how to attend is on the LAMMPS homepage at
[Public Meetings](https://www.lammps.org/community/meetings/).

### 2.3 LAMMPS Slack

The LAMMPS project uses slack for day-to-day developer communication at
[lammps.slack.com](https://lammps.slack.com).  Invitations to the LAMMPS
slack instance can be requested by sending an email to slack@lammps.org.

The LAMMPS Slack is low volume and intended for developers only.

### 2.4 LAMMPS Github

LAMMPS [issues](https://github.com/lammps/lammps/issues) and [pull
requests](https://github.com/lammps/lammps/pulls), are tracked on the
[LAMMPS GitHub page](https://github.com/lammps/lammps).

### 2.5 Developer Email Address Alias

Anybody without access to the websites above or who wants to contact the
LAMMPS developers privately, can send an email to developers@lammps.org
and it will be forwarded to the current LAMMPS core developers.

---

## 3. Announcements

Announcements will be posted on the LAMMPS website
[lammps.org](https://www.lammps.org/) and on
[MatSci.org/lammps](https://matsci.org/lammps).

---

## 4. Teams within the LAMMPS project

<br>

### 4.1 Administrative Access to the Github Repository

Only members of the TSC may have administrative access to the LAMMPS
repository.  TSC members will be given access individually as needed.
There must be *at least* two TSC members with the GitHub "Owner"
permission at all times.  This list is not generally published for
operational security reasons.

### 4.2 Other developer groups within the LAMMPS project

There are currently three more categories of GitHub users with
additional status or permissions in the LAMMPS GitHub project.

- **Contributors** are GitHub users contributing non-trivial changes to
  LAMMPS.  They are sent an invitation to obtain a "Collaborator" status
  on the project by GitHub and may be given "Triage" permission on
  issues and pull requests.

- **Core Developers** have write access to public and private
  repositories in the project except for protected branches.  Two
  approvals from *Core developers* are required to merge a pull request;
  the step of merging by a *Maintainer* counts as an implicit approval.
  *Core Developers* may also request changes to pull requests which will
  block merging them until they approve.

- **Maintainers** are *Core Developers* authorized to merge pull
  requests and make administrative changes to the GitHub project and any
  servers operated on behalf of the project.

Final decisions on membership in any of these teams are made by the TSC.
*All* participants in the LAMMPS project are collectively referred to as
**Collaborators**.

---

A [copy of this document](https://github.com/lammps/lammps/blob/develop/.github/GOVERNANCE.md)
is also available in the LAMMPS git repository as `.github/GOVERNANCE.md`.
