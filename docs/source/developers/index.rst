.. _developers:

====================================
Developers: installation from source
====================================

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :caption: Contents:
   
   pyfabm-dev
   fabm0d-dev
   gotm-dev


Obtaining source code
~~~~~~~~~~~~~~~~~~~~~

To get the ERSEM and FABM source codes:

.. code-block:: bash

   git clone https://github.com/pmlmodelling/ersem.git
   git clone https://github.com/fabm-model/fabm.git

This locally creates ``ersem`` and ``fabm`` directories with source
code.

For water column models, you will need the GOTM source code:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/gotm-model/code.git gotm

To build the code, you will need: \* `a recent Fortran
compiler <https://github.com/fabm-model/fabm/wiki/Building-and-installing#supported-compilers>`__
\* `cmake <https://www.cmake.org>`__ 3.0 or higher. First check whether
you have that installed: run ``cmake --version`` on the command line.

.. note::
    It is not necessary to use ``pip`` or ``apt`` to install dependences. We 
    suggest that `conda <https://docs.conda.io/en/latest/>`__ or 
    `brew <https://brew.sh/>`__ would also work well.

