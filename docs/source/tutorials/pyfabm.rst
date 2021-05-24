.. _pyfabm:

#############################
Python front-end (``pyfabm``)
#############################

``FABM`` comes with a ``pyfabm`` package for the ``python`` programming language.
This package enables you to access 
`FABM's <https://github.com/fabm-model/fabm/wiki>`__ biogeochemical models directly 
from ``python``. For instance, you can enumerate a model's parameters and variables, 
or obtain temporal derivatives and diagnostics for any given environment and 
model state. In combination with a ``python``-based time integration scheme 
(e.g., ``scipy.integrate.odeint``), this also allows you to perform model 
simulations. More information about ``pyfabm`` can be found on the 
`FABM wiki page <https://github.com/fabm-model/fabm/wiki/python>`__. 
Below you can find brief instructions on how to use ``pyfabm`` with ERSEM.

Running ``pyfabm`` tutorial
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To demonstrate how ``pyfabm`` can be used with ERSEM, we will calculate the 
oxygen saturation concentration as a function of temperature and salinity 
using the method implemented in ERSEM, which is taken from :cite:`WEISS1970`.

.. note::
    The following script requires ``matplotlib`` to be installed. This can
    easily be done via ``pip`` in the following way:

    .. code-block:: bash
        
        python -m pip install matplotlib

With ``pyfabm`` installed, this can be achieved using the following code:

.. literalinclude:: ../../../github-actions/pyfabm-ersem/pyfabm_tut.py
    :language: python
    :linenos:

.. image:: ../../images/pyfabm_fig.png
   :alt: Example ``pyfabm`` output
   :width: 100.0%

In additions, example Jupyter notebooks that use the Python front-end can 
be found in ``<FABMDIR>/testcases/python`` and more examples can be found 
on `the FABM wiki <https://github.com/fabm-model/fabm/wiki/python>`__.

