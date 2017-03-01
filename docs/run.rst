.. _run:

Usage
====================================

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

.. argparse::
    :module: lib.main
    :func: define_parser
    :prog: pvacseq run
