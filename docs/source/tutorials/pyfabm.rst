.. _pyfabm:

`pyfabm`
~~~~~~~~

To demonstrate how `pyfabm` can be used with ERSEM, we will calculate the 
oxygen saturation concentration as a function of temperature and salinity 
using the method implemented in ERSEM, which is taken from Weiss (1970).
ith `pyfabm` installed, this can be achieved using the following code:

.. literalinclude:: ../../../github-actions/pyfabm-ersem/pyfabm-tut.py
    :language: python
    :linenos:

.. image:: ../../images/pyfabm_fig.png
   :alt: Example `pyfabm` output
   :width: 100.0%

In additions, example Jupyter notebooks that use the Python front-end can 
be found in ``<FABMDIR>/testcases/python`` and more examples can be found 
on `the FABM wiki <https://github.com/fabm-model/fabm/wiki/python>`__.

