---
title: "Charter"
description: "Scope and goals of the LAMMPS project"
---

<sub>(this page and its content is currently a work in progress and not final)</sub>

The "Charter" for the LAMMPS Project, part of a Series of LF Projects

Adopted August xx, 2026.

This Charter sets forth the responsibilities and procedures for
technical contribution to, and oversight of the LAMMPS open source
project, which has been established as LAMMPS Project a Series of LF
Projects (the "Project").  LF Projects, LLC ("LF Projects") is a
Delaware series limited liability company.  All contributors (including
Maintainers, Core Developers, and other technical positions) and other
participants in the Project (collectively, "Collaborators") must comply
with the terms of this Charter.

## 1. Mission and Scope of the Project

The mission of the Project is to:

- Develop and support the LAMMPS code that is used for particle
  simulations with a focus on large-scale classical molecular
  dynamics for materials running on multi-node parallel computing
  platforms ("HPC Clusters");

- Support a diverse community of users in materials science, materials
   chemistry, condensed matter physics, and related fields;

- Balance the benefits of modern software engineering practices
  against the need for the code to be understandable and modifiable by
  domain scientists with limited software development experience;

- Continue to adapt the code base to support both:

  - new HPC Cluster hardware architectures, including conventional
     CPU and GPU architectures, as well emerging architectures; and

  - new physics models, including those built on emerging machine
    learning and artificial intelligence technologies;

- Promote the co-design of physics models and algorithms to drive the
  development of Pareto-optimal implementations that balance physical
  realism against computational cost;

- Maintain a comprehensive, secure software distribution, and

- Provide documentation, training, and user support.

The scope of the Project includes collaborative development under the
Project License (as defined herein) supporting the mission, including
documentation, testing, integration and the creation of other artifacts
that aid the development, deployment, operation or adoption of the open
source project.

## 2. Technical Steering Committee

The Technical Steering Committee (the "TSC") will be responsible for all
technical oversight of the open source Project.

- The TSC voting members are initially the people set forth within the
  [Governance](https://www.lammps.org/about/governance) page.  The TSC may choose an
  alternative approach for determining the voting members of the TSC,
  and any such alternative approach will be documented in the
  [Governance](https://www.lammps.org/about/governance) page.  Any meetings of the Technical
  Steering Committee are intended to be open to the public, and can be
  conducted electronically, via teleconference, or in person.

- TSC projects generally will involve Contributors, Core Developers, and
  Maintainers.  The TSC may adopt or modify roles so long as the roles
  are documented in the [Governance](https://www.lammps.org/about/governance) page. Unless
  otherwise documented roles are:

  - **Contributors** include anyone in the technical community that
    contributes code, documentation, or other technical artifacts to
    the Project;

  - **Core Developers** are Contributors who have contributed code,
    documentation, and other technical artifacts over a multi-year
    period and/or have unique expertise in specific aspects of the
    Project;

  - **Maintainers** are Core Developers who have earned the ability to
    merge changes to ("commit") source code, documentation or other
    technical artifacts in a project’s repository;

  - Core Developers are **recruited** by the existing Core Developers; and

  - Maintainers are **appointed** by the TSC. Anyone, but typically a
    Core Developer, may become a Maintainer by a majority approval of
    the existing TSC. A Maintainer may be removed by a majority
    approval of the TSC.

- Participation in the Project through becoming a Contributor, Core
  Developer, or Maintainer is open to anyone so long as they abide by
  the terms of this Charter.

- The TSC may (1) establish workflow procedures for the submission,
  approval, and closure/archiving of projects, (2) set requirements for
  the promotion of Core Developers to Maintainer status, as applicable,
  and (3) amend, adjust, refine and/or eliminate the roles of
  Contributors, Core Developers, and Maintainers, and create new roles,
  and publicly document any TSC roles, as it sees fit.

- The TSC may elect a TSC Chair, who will preside over meetings of the
  TSC and will serve until their resignation or replacement by the TSC.

- Responsibilities: The TSC will be responsible for all aspects of
  oversight relating to the Project, which may include:

  1. coordinating the technical direction of the Project;

  2. approving project or system proposals (including, but not limited
     to, incubation, deprecation, and changes to a sub-project’s scope);

  3. establishing community norms, workflows, issuing releases, and
     security issue reporting policies;

  4. approving and implementing policies and processes for contributing
    (published in a
    [.github/CONTRIBUTING.md](https://github.com/lammps/lammps/blob/develop/.github/CONTRIBUTING.md)
    file in the LAMMPS GitHub repository) and coordinating with the
    series manager of the Project (as provided for in the Series
    Agreement, the "Series Manager") to resolve matters or concerns that
    may arise as set forth in Section 7 of this Charter;

  5. discussions, seeking consensus, and where necessary, voting on
     technical matters relating to the code base; and

  6. coordinating any marketing, events, or communications regarding the Project.

  7. All decisions are made by seeking broad consensus among the TSC
     members.  This happens primarily through:

       a) discussion and review on GitHub pull requests and issues;

       b) monthly developer meetings;

       c) quarterly TSC meetings

     If consensus cannot be reached through discussion, the decision
     will be put to a vote at a TSC meeting according to section 3.

## 3. TSC Voting

- While the Project aims to operate as a consensus-based community, if
  any TSC decision requires a vote to move the Project forward, the
  voting members of the TSC will vote on a one vote per voting member
  basis.

- Quorum for TSC meetings requires at least fifty percent (50\%) of all
  voting members of the TSC to be present. The TSC may continue to meet
  if quorum is not met but will be prevented from making any decisions
  at the meeting.

- Decisions by vote at a meeting require a majority vote of TSC members
  in attendance, provided quorum is met. Decisions made by electronic
  vote without a meeting require a majority vote of all voting members
  of the TSC.

- In the event a vote cannot be resolved by the TSC, any voting member
  of the TSC may refer the matter to the Series Manager for assistance
  in reaching a resolution.

## 4. Compliance with Policies

- This Charter is subject to the Series Agreement for the Project and
  the Operating Agreement of LF Projects. Contributors will comply with
  the policies of LF Projects as may be adopted and amended by LF
  Projects, including, without limitation the policies listed at
  [https://lfprojects.org/policies/](https://lfprojects.org/policies/)

- The TSC may adopt a code of conduct ("CoC") for the Project, which is
  subject to approval by the Series Manager. In the event that a
  Project-specific CoC has not been approved, the LF Projects Code of
  Conduct listed at
  [https://lfprojects.org/policies](https://lfprojects.org/policies)
  will apply for all Collaborators in the Project.

- When amending or adopting any policy applicable to the Project, LF
  Projects will publish such policy, as to be amended or adopted, on its
  web site at least 30 days prior to such policy taking effect;
  provided, however, that in the case of any amendment of the Trademark
  Policy or Terms of Use of LF Projects, any such amendment is effective
  upon publication on LF Project’s web site.

- All Collaborators must allow open participation from any individual or
  organization meeting the requirements for contributing under this
  Charter and any policies adopted for all Collaborators by the TSC,
  regardless of competitive interests. Put another way, the Project
  community must not seek to exclude any participant based on any
  criterion, requirement, or reason other than those that are reasonable
  and applied on a non-discriminatory basis to all Collaborators in the
  Project community.

- The Project will operate in a transparent, open, collaborative, and ethical manner at
  all times. The output of all Project discussions, proposals, timelines, decisions,
  and status should be made open and easily visible to all. Any potential violations of
  this requirement should be reported immediately to the Series Manager.

## 5. Community Assets

- LF Projects will hold title to all trade or service marks used by the Project
  ("Project Trademarks"), whether based on common law or registered rights. Project
  Trademarks will be transferred and assigned to LF Projects to hold on behalf of the
  Project. Any use of any Project Trademarks by Collaborators in the Project will be in
  accordance with the license from LF Projects and inure to the benefit of LF Projects.

- The Project will, as permitted and in accordance with such license from LF Projects,
  develop and own all Project GitHub and social media accounts, and domain name
  registrations created by the Project community.

- Under no circumstances will LF Projects be expected or required to undertake any
  action on behalf of the Project that is inconsistent with the tax-exempt status or
  purpose, as applicable, of the Joint Development Foundation or LF Projects, LLC.

## 6. General Rules and Operations

The Project will:

  1. engage in the work of the Project in a professional manner consistent with
     maintaining a cohesive community, while also maintaining the goodwill and esteem
     of LF Projects, Joint Development Foundation and other partner organizations in
     the open source community; and

  2. respect the rights of all trademark owners, including any branding and trademark
     usage guidelines.

## 7. Intellectual Property Policy

- Collaborators acknowledge that the copyright in all new contributions will be
  retained by the copyright holder as independent works of authorship and that no
  contributor or copyright holder will be required to assign copyrights to the Project.

- All contributions to the Project are subject to the following:

    1. All new inbound code contributions to the Project must be made using both
       (1) the GNU General Public License version 2
       [(SPDX: GPL-2.0-only)](https://spdx.org/licenses/GPL-2.0-only.html)
       and (2) GNU Lesser General Public License version 2.1
       [(SPDX: LGPL-2.1-only)](https://spdx.org/licenses/LGPL-2.1-only.html)
       (together, the "Project License").

    2. All new inbound non-trivial code contributions must also be
       accompanied by the filled out pull request template where the
       authors are asked to:
         - identify themselves and their affiliation(s)
         - provide a (persistent) email address for future communication
           with the main author(s)
         - agree to the licensing terms for the contributed code
         - disclose usage of AI coding tools

       Pull requests without the filled out template, with the exception
       of trivial changes or pull requests required for project management
       (e.g. to set version strings for a release), will be **rejected**.

    3. All outbound code will be made available under the Project License.

    4. Documentation will be received and made available by the Project under the
      [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

- Contributed source files should contain license information indicating the open source
  license or licenses pertaining to the file.

## 8. Amendments

This charter may be amended by a two-thirds vote of the entire TSC and
is subject to approval by LF Projects.

-----

A [copy of this document](https://github.com/lammps/lammps/blob/develop/.github/CHARTER.md)
is also available in the LAMMPS git repository as `.github/CHARTER.md`.
