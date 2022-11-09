.. _gotm-dev:

~~~~~~~~~~~~~~
``gotm-ersem``
~~~~~~~~~~~~~~


To run the install script, you will need to have ``netCDF`` installed.
An example of how to do this is here:

.. literalinclude:: ../../../github-actions/gotm-fabm-ersem/gotm-fabm-ersem-dep-debian.sh
    :language: bash
    :linenos:

To install GOTM-FABM-ERSEM, we suggest you use the following script below

.. literalinclude:: ../../../github-actions/gotm-fabm-ersem/gotm-fabm-ersem-build.sh
    :language: bash
    :linenos:


If you experience NetCDF issues when running ``make install``, see :ref:`trouble`.

Now you should have a GOTM executable with FABM and ERSEM support at
``~/local/gotm/bin/gotm``.
