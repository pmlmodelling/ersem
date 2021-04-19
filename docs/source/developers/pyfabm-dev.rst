.. _pyfabm-dev:

~~~~~~~~~~~~~~~~
``pyfabm-ersem``
~~~~~~~~~~~~~~~~

To build the ``python`` driver and use it with ``ERSEM``, you must first obtain 
copies of the ``FABM`` and ``ERSEM`` codes. The latest version of each can be 
obtained by cloning their respective code repositories using 
`git <https://git-scm.com/>`__. The instructions of how to do this are found
:ref:`here <developers>`.

To run the install script you will need to have ``wheel`` and ``numpy`` installed.
This can be done via:

.. literalinclude:: ../../../github-actions/pyfabm-ersem/pyfabm-ersem-dep-debian.sh
    :language: bash
    :linenos:

To install PyFABM-ERSEM, we suggest you use the following script below

.. literalinclude:: ../../../github-actions/pyfabm-ersem/pyfabm-ersem-build.sh
    :language: bash
    :linenos:
