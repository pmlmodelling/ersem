.. _ersemworkflow:

##############
ERSEM Workflow
##############

As a rule of thumb, ERSEM developments should be made in the `public ERSEM <https://github.com/pmlmodelling/ersem>`_
with the exception of specific project developments.

Issue tracking
--------------

All developments, either private or public, are required to create an 
`issue/issues <https://github.com/pmlmodelling/ersem/issues>`_ describing the problem or
development required.

For larger project developments that require multiple issues, it is often preferable to create a
`milestone <https://github.com/pmlmodelling/ersem/milestones>`_ which you can then assign multiple
issues too.

Private developments
--------------------

The only private developments would come from 
specific projects - these require significant scientific changes to the model,
and thus, are isolated from the public version of the code.

Workflow for private developments
+++++++++++++++++++++++++++++++++

The workflow for private developments is as follows:

1. Create an issue or issues in the `public ERSEM <https://github.com/pmlmodelling/ersem>`_ 
   repo based on the development required.
2. Create a `branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`_
   from the `private ERSEM <https://github.com/pmlmodelling/ersem-dev>`_ repo. The naming
   convention of the branch is as follows `XX-[description]` where `XX` is the issue ID (this can be found
   after creation of the issue) and `[description]` is a one- or two- word description to identify the branch
   with the corresponding issue.
3. Commit all changes outlined in the issue to this branch and push branch to the private ERSEM repo.
4. Create `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_
   to merge the branch you have made the changes on into the `master` branch. In the creation of the pull request you
   will need to assign specific people to review the changes you have made before it can be merged into the master branch.
5. Once the updates have been approved by the reviewer/reviewers they will merge the development into the master branch.
6. When the private master branch is updated you will need to talk to the ERSEM repo manager to merge those changes back
   to the public repo.

Public developments
-------------------

The following come under public developments:

1. Documentation - new descriptions and updates to current documentation, these include new ERSEM tutorials or 
   fixes to current docs
2. Bug fixes - these require immediate attention when identified within the code base
3. Public project developments - if possible, small developments and additions to the code shall be added publicly,
   these could include different options or parameterisations of the model.

Workflow for public developments
++++++++++++++++++++++++++++++++

For public developments, the workflow is pretty simple:

1. Create an issue or issues in the `public ERSEM <https://github.com/pmlmodelling/ersem>`_ 
   repo based on the development required.
2. Create a `branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`_
   from the `public ERSEM <https://github.com/pmlmodelling/ersem>`_ repo. The naming
   convention of the branch is as follows `XX-[description]` where `XX` is the issue ID (this can be found
   after creation of the issue) and `[description]` is a one- or two- word description to identify the branch
   with the corresponding issue.
3. Commit all changes outlined in the issue to this branch and push branch to the ERSEM repo.
4. Create `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_
   to merge the branch you have made the changes on into the `master` branch. In the creation of the pull request you
   will need to assign specific people to review the changes you have made before it can be merged into the master branch.
5. Once the updates have been approved by the reviewer/reviewers they will merge the development into the master branch.
6. When the public master branch is updated you will need to talk to the ERSEM repo manager to merge those changes back
   to the private repo.

Notes for reviewers
-------------------

When `reviewing pull requests <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews>`_ there are a few things to check:

1. Ensure that the new changes do not adversely affect ERSEM. **Note** some modifications will change the results
   produced by the model, however, it is the author of the changes and reviewer's responsibility to decide
   whether these changes go into ERSEM. If unsure, please contact the ERSEM repository manager.
2. Ensure that the changes made to the code correctly address the issue created.
3. For the public ERSEM, make sure that the tests continue to pass.

.. note::

    Any questions about the workflow, please contact the ERSEM repository manager.

