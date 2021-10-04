.. _tutorials:


ERSEM tutorials
===============

ERSEM is mainly run with 3D hydrodynamic models such as 
`NEMO <https://www.nemo-ocean.eu/>`__  and 
`FVCOM <http://fvcom.smast.umassd.edu/fvcom/>`__. However, due to its
flexiblilty, it can also be run as a 0D box model or within a 1D water column,
via `GOTM <https://gotm.net/portfolio/>`__. Parts of the library can also 
be called directly from Python through the ``pyfabm`` interface.
Here we describe several tutorials for the common use cases for ERSEM.

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :caption: Contents:

   pyfabm
   fabm0d
   gotm
   fvcom
   nemo

.. _conda-install:

``conda`` installation
~~~~~~~~~~~~~~~~~~~~~~

We recommend using ``conda`` to obtain the latest version of ``ersem``. The
conda package includes interfaces for :ref:`pyfabm`, :ref:`fabm0d` and 
:ref:`gotm` tutorials, and can be installed as follows:

.. code-block:: bash

    conda create -n ersem-tut -y
    conda activate ersem-tut
    conda install -c pmlmodelling ersem -y


To install ``ersem`` from source please look at the 
:ref:`developers documentation <developers>`.
