.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

.. _pvacvector_run:

Usage
====================================

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

.. topic:: For usage instructions run
   
   ``pvavector run --help``

.. .. argparse::
        :module: tools.pvacvector.run
        :func: define_parser
        :prog: pvacvector run
