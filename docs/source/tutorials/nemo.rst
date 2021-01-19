.. _nemo:


NEMO + FABM + ERSEM
~~~~~~~~~~~~~~~~~~~

This requires a customised NEMO codebase with the FABM coupler
integrated. To obtain access to this code, please
`register <https://pml.ac.uk/Modelling_at_PML/Access_Code>`__.

FABM-ERSEM needs to be compiled separately before the compilation of
NEMO. Usually, the following suffices for this:

::

   mkdir -p ~/build/nemo
   cd ~/build/nemo
   cmake <SOURCEDIR>/fabm -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<SOURCEDIR>/ersem
   make install

This will create the library in the standard folder
``~/local/fabm/nemo/lib`` where NEMO-FABM will look for linking to NEMO.
Module files (.mod) are placed at ``~/local/fabm/nemo/include``.


