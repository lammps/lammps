# Contributing to LAMMPS via GitHub

Thank your for considering to contribute to the LAMMPS software project.

The following is a set of guidelines as well as explanations of policies and work flows for contributing to the LAMMPS molecular dynamics software project. These guidelines focus on submitting issues or pull requests on the LAMMPS GitHub project.

Thus please also have a look at:
* [The guide for submitting new features in the LAMMPS manual](https://www.lammps.org/doc/Modify_contribute.html)
* [The guide on programming style and requirement in the LAMMPS manual](https://www.lammps.org/doc/Modify_style.html)
* [The GitHub tutorial in the LAMMPS manual](http://lammps.sandia.gov/doc/Howto_github.html)

## Table of Contents

[I don't want to read this whole thing, I just have a question!](#i-dont-want-to-read-this-whole-thing-i-just-have-a-question)

[How Can I Contribute?](#how-can-i-contribute)
* [Discussing How To Use LAMMPS](#discussing-how-to-use-lammps)
* [Reporting Bugs](#reporting-bugs)
* [Suggesting Enhancements](#suggesting-enhancements)
* [Contributing Code](#contributing-code)

[GitHub Work flows](#github-workflows)
* [Issues](#issues)
* [Pull Requests](#pull-requests)

__

## I don't want to read this whole thing I just have a question!

> **Note:** Please do not file an issue to ask a general question about LAMMPS, its features, how to use specific commands, or how perform simulations or analysis in LAMMPS. Instead post your question to either the ['lammps-users' mailing list](https://lammps.sandia.gov/mail.html) or the [LAMMPS Material Science Discourse forum](https://matsci.org/lammps). You do not need to be subscribed to post to the list (but a mailing list subscription avoids having your post delayed until it is approved by a mailing list moderator). Most posts to the mailing list receive a response within less than 24 hours. Before posting to the mailing list, please read the [mailing list guidelines](https://lammps.sandia.gov/guidelines.html). Following those guidelines will help greatly to get a helpful response. Always mention which LAMMPS version you are using. The LAMMPS forum was recently created as part of a larger effort to build a materials science community and have discussions not just about using LAMMPS. Thus the forum may be also used for discussions that would be off-topic for the mailing list. Those will just have to be posted to a more general category.

## How Can I Contribute?

There are several ways how you can actively contribute to the LAMMPS project: you can discuss compiling and using LAMMPS, and solving LAMMPS related problems with other LAMMPS users on the lammps-users mailing list or the forum, you can report bugs or suggest enhancements by creating issues on GitHub (or posting them to the lammps-users mailing list or posting in the LAMMPS Materials Science Discourse forum), and you can contribute by submitting pull requests on GitHub or e-mail your code
to one of the [LAMMPS core developers](https://lammps.sandia.gov/authors.html). As you may see from the aforementioned developer page, the LAMMPS software package includes the efforts of a very large number of contributors beyond the principal authors and maintainers.

### Discussing How To Use LAMMPS

The LAMMPS mailing list is hosted at SourceForge. The mailing list began in 2005, and now includes tens of thousands of messages in thousands of threads. LAMMPS developers try to respond to posted questions in a timely manner, but there are no guarantees. Please consider that people live in different timezone and may not have time to answer e-mails outside of their work hours.
You can post to list by sending your email to lammps-users at lists.sourceforge.net (no subscription required), but before posting, please read the [mailing list guidelines](https://lammps.sandia.gov/guidelines.html) to maximize your chances to receive a helpful response.

Anyone can browse/search previous questions/answers in the archives. You do not have to subscribe to the list to post questions, receive answers (to your questions), or browse/search the archives. You **do** need to subscribe to the list if you want emails for **all** the posts (as individual messages or in digest form), or to answer questions yourself. Feel free to sign up and help us out! Answering questions from fellow LAMMPS users is a great way to pay back the community for providing you a useful tool for free, and to pass on the advice you have received yourself to others. It improves your karma and helps you understand your own research better.

If you post a message and you are a subscriber, your message will appear immediately. If you are not a subscriber, your message will be moderated, which typically takes one business day. Either way, when someone replies the reply will usually be sent to both, your personal email address and the mailing list. When replying to people, that responded to your post to the list, please always included the mailing list in your replies (i.e. use "Reply All" and **not** "Reply"). Responses will appear on the list in a few minutes, but it can take a few hours for postings and replies to show up in the SourceForge archive. Sending replies also to the mailing list is important, so that responses are archived and people with a similar issue can search for possible solutions in the mailing list archive.

The LAMMPS Materials Science Discourse forum was created recently to facilitate discussion not just about LAMMPS and as part of a larger effort  towards building a materials science community. The forum contains a read-only sub-category with the continually updated mailing list archive, so you won't miss anything by joining only the forum and not the mailing list.

### Reporting Bugs

While developers writing code for LAMMPS are careful to test their code, LAMMPS is such a large and complex software, that it is impossible to test for all combinations of features under all normal and not so normal circumstances. Thus bugs do happen, and if you suspect, that you have encountered one, please try to document it and report it as an [Issue](https://github.com/lammps/lammps/issues) on the LAMMPS GitHub project web page. However, before reporting a bug, you need to check whether this is something that may have already been corrected. The [Latest Features and Bug Fixes in LAMMPS](https://lammps.sandia.gov/bug.html) web page lists all significant changes to LAMMPS over the years. It also tells you what the current latest development version of LAMMPS is, and you should test whether your issue still applies to that version.

When you click on the green "New Issue" button, you will be provided with a text field, where you can enter your message. That text field with contain a template with several headlines and some descriptions. Keep the headlines that are relevant to your reported potential bug and replace the descriptions with the information as suggested by the descriptions.
You can also attach small text files (please add the file name extension `.txt` or it will be rejected), images, or small compressed text files (using gzip, do not use RAR or 7-ZIP or similar tools that are uncommon outside of Windows machines). In many cases, bugs are best illustrated by providing a small input deck (do **not** attach your entire production input, but remove everything that is not required to reproduce the issue, and scale down your system size, that the resulting calculation runs fast and can be run on small desktop quickly).

To be able to submit an issue on GitHub, you have to register for an account (for GitHub in general). If you do not want to do that, or have other reservations against submitting an issue there, you can - as an alternative and in decreasing preference - either send an e-mail to the lammps-users mailing list, the original authors of the feature that you suspect to be affected, or one or more of the core LAMMPS developers.

### Suggesting Enhancements

The LAMMPS developers welcome suggestions for enhancements or new features. These should be submitted using the [GitHub Issue Tracker](https://github.com/lammps/lammps/issues) of the LAMMPS project. This is particularly recommended, when you plan to implement the feature or enhancement yourself, as this allows to coordinate in case there are other similar or conflicting ongoing developments.
The LAMMPS developers will review your submission and consider implementing it. Whether this will actually happen depends on many factors: how difficult it would be, how much effort it would take, how many users would benefit from it, how well the individual developer would understand the underlying physics of the feature, and whether this is a feature that would fit into a software like LAMMPS, or would be better implemented as a separate tool. Because of these factors, it matters how well the suggested enhancement is formulated and the overall benefit is argued convincingly.

To be able to submit an issue on GitHub, you have to register for an account (for GitHub in general). If you do not want to do that, or have other reservations against submitting an issue there, you can - as an alternative - send an e-mail to the lammps-users mailing list.

### Contributing Code

We encourage users to submit new features or modifications for LAMMPS. Instructions, guidelines, requirements,
and recommendations are in the following sections of the LAMMPS manual:
* [The guide for submitting new features in the LAMMPS manual](https://lammps.sandia.gov/doc/Modify_contribute.html)
* [The guide on programming style and requirement in the LAMMPS manual](https://lammps.sandia.gov/doc/Modify_contribute.html)
* [The GitHub tutorial in the LAMMPS manual](http://lammps.sandia.gov/doc/Howto_github.html)


## GitHub Workflows

This section briefly summarizes the steps that will happen **after** you have submitted either an issue or a pull request on the LAMMPS GitHub project page. 

### Issues

After submitting an issue, one or more of the LAMMPS developers will review it and categorize it by assigning labels. Confirmed bug reports will be labeled `bug`; if the bug report also contains a suggestion for how to fix it, it will be labeled `bugfix`; if the issue is a feature request, it will be labeled `enhancement`. Other labels may be attached as well, depending on which parts of the LAMMPS code are affected. If the assessment is, that the issue does not warrant any changes, the `wontfix` label will be applied and if the submission is incorrect or something that should not be submitted as an issue, the `invalid` label will be applied. In both of the last two cases, the issue will then be closed without further action.

For feature requests, what happens next is that developers may comment on the viability or relevance of the request, discuss and make suggestions for how to implement it. If a LAMMPS developer or user is planning to implement the feature, the issue will be assigned to that developer. For developers, that are not yet listed as LAMMPS project collaborators, they will receive an invitation to be added to the LAMMPS project as a collaborator so they can get assigned. If the requested feature or enhancement is implemented, it will be submitted as a pull request, which will contain a reference to the issue number. And once the pull request is reviewed and accepted for inclusion into LAMMPS, the issue will be closed. For details on how pull requests are processed, please see below. Feature requests may be labeled with `volunteer_needed` if none of the LAMMPS developers has the time and the required knowledge implement the feature.

For bug reports, the next step is that one of the core LAMMPS developers will self-assign to the issue and try to confirm the bug. If confirmed, the `bug` label and potentially other labels are added to classify the issue and its impact to LAMMPS. Otherwise the `unconfirmed` label will be applied and some comment about what was tried to confirm the bug added. Before confirming, further questions may be asked or requests for providing additional input files or details about the steps required to reproduce the issue. Any bugfix will be submitted as a pull request (more about that below) and since most bugs require only local changes, the bugfix may be included in a pull request specifically set up to collect such local bugfixes or small enhancements. Once the bugfix is included in the master branch, the issue will be closed.

### Pull Requests

Pull requests are the **only** way that changes get made to the LAMMPS distribution.  So also the LAMMPS core developers will submit pull requests for their own changes and discuss them on GitHub.  Thus if you submit a pull request it will be treated in a similar fashion. When you submit a pull request you may opt to submit a "Draft" pull request.  That means your changes are visible and will be subject to testing, but reviewers will not be (auto-)assigned and comments will take into account that this is not complete. On the other hand, this is a perfect way to ask the LAMMPS developers for comments on non-obvious changes and get feedback and possible suggestions for improvements or recommendations about what to avoid.
Immediately after the submission, the LAMMPS continuing integration server at ci.lammps.org will download your submitted branch and perform a number of tests: it will tests whether it compiles cleanly under various conditions, it will also do a check on whether your included documentation translates cleanly and run some unit tests and other checks. Whether these tests are successful or fail will be recorded.  If a test fails, please inspect the corresponding output on the CI server and take the necessary steps, if needed, so that the code can compile cleanly again.  The test will be re-run each time the pull request is updated with a push to the remote branch on GitHub.  If you are unsure about what you need to change, ask a question in the discussion area of the pull request.
Next a LAMMPS core developer will self-assign and do an overall technical assessment of the submission.  If you submitted a draft pull request, this will not happen unless you mark it "ready for review".  If you are not yet invited as a LAMMPS collaborator, and your contribution seems significant, you may also receive an invitation for collaboration on the LAMMPS repository.  As part of the assessment, the pull request will be categorized with labels. There are two special labels: `needs_work` (indicates that work from the submitter of the pull request is needed) and `work_in_progress` (indicates, that the assigned LAMMPS developer will make changes, if not done by the contributor who made the submit). 
You may also receive comments and suggestions on the overall submission or specific details and on occasion specific requests for changes as part of the review. If permitted, also additional changes may be pushed into your pull request branch or a pull request may be filed in your LAMMPS fork on GitHub to include those changes.
The LAMMPS developer may then decide to assign the pull request to another developer (e.g. when that developer is more knowledgeable about the submitted feature or enhancement or has written the modified code). It may also happen, that additional developers are requested to provide a review and approve the changes. For submissions, that may change the general behavior of LAMMPS, or where a possibility of unwanted side effects exists, additional tests may be requested by the assigned developer.
If the assigned developer is satisfied and considers the submission ready for inclusion into LAMMPS, the pull request will receive approvals and be merged into the master branch by one of the core LAMMPS developers. After the pull request is merged, you may delete the feature branch used for the pull request in your personal LAMMPS fork.  The minimum requirement to merge a pull request is that all automated tests have to pass and at least one LAMMPS developer has approved integrating the submitted code. Since the approver will not be the person merging a pull request, you will have at least two LAMMPS developers that looked at your contribution.
Since the learning curve for git is quite steep for efficiently managing remote repositories, local and remote branches, pull requests and more, do not hesitate to ask questions, if you are not sure about how to do certain steps that are asked of you. Even if the changes asked of you do not make sense to you, they may be important for the LAMMPS developers. Please also note, that these all are guidelines and nothing set in stone. So depending on the nature of the contribution, the work flow may be adjusted.

