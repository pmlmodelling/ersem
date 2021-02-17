.. _tutorials:


ERSEM use cases
===============

In this section, we give several tutorials of common use cases for ERSEM.

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :caption: Contents:

   pyfabm
   fabm0d
   gotm
   nemo

.. note::
    The tutorials presented here are run via 
    `github actions <https://github.com/pmlmodelling/ersem/actions>`__ 
    on every pull request.

Obtaining source code
~~~~~~~~~~~~~~~~~~~~~

Now get the ERSEM, FABM and GOTM source codes:

.. code-block:: bash

   git clone https://github.com/pmlmodelling/ersem.git
   git clone https://github.com/fabm-model/fabm.git
   git clone --recurse-submodules https://github.com/gotm-model/code.git gotm

This locally creates ``gotm``, ``fabm``, ``ersem`` directories with source
code.

Note: the above gets you ERSEM’s public stable release. Developers can
check out the the developers’ version by substituting
``git@gitlab.ecosystem-modelling.pml.ac.uk:edge/ersem.git`` for
``https://github.com/pmlmodelling/ersem.git``.

To build the code, you will need: \* `a recent Fortran
compiler <https://github.com/fabm-model/fabm/wiki/Building-and-installing#supported-compilers>`__
\* `cmake <https://www.cmake.org>`__ 3.0 or higher. First check whether
you have that installed: run ``cmake --version`` on the command line.

.. note::
    It is not necessary to use ``pip`` or ``apt`` to install dependences. We 
    suggest that `conda <https://docs.conda.io/en/latest/>`__ or 
    `brew <https://brew.sh/>`__ would also work well.


