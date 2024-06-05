.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Additional Commands
===================

To make using pVACfuse easier, several convenience methods are included in the package.

.. _pvacfuse_example_data:

Download Example Data
---------------------

.. program-output:: pvacfuse download_example_data -h

.. .. argparse::
    :module: lib.download_example_data
    :func: define_parser
    :prog: pvacfuse download_example_data

List Valid Alleles
------------------

.. program-output:: pvacfuse valid_alleles -h

.. .. argparse::
    :module: lib.valid_alleles
    :func: define_parser
    :prog: pvacfuse valid_alleles

.. _valid_algorithms:

List Valid Algorithms
---------------------

.. program-output:: pvacfuse valid_algorithms -h

.. .. argparse::
    :module: lib.valid_algorithms
    :func: define_parser
    :prog: pvacfuse valid_algorithms

List Allele-Specific Cutoffs
----------------------------

.. program-output:: pvacfuse allele_specific_cutoffs -h
