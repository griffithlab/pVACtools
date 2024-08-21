.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Additional Commands
===================

To make using pVACvector easier, several convenience methods are included in the package.

Creating Vector Visualization
-----------------------------

By default, pVACvector will create a visualization of the vector design
result. For this to happen, the DISPLAY environment variable has to be set.
This is often not the case, for example, when running pVACvector on a compute
cluster. We provide this convenience method to create the visualization
graphic from a successful pVACvector result FASTA file on any machine that has
the DISPLAY environment variable set.

.. program-output:: pvacvector visualize -h

.. _pvacvector_example_data:

Download Example Data
---------------------

.. program-output:: pvacvector download_example_data -h

.. .. argparse::
    :module: lib.download_example_data
    :func: define_parser
    :prog: pvacfuse download_example_data

List Valid Alleles
------------------

.. program-output:: pvacvector valid_alleles -h

.. .. argparse::
    :module: lib.valid_alleles
    :func: define_parser
    :prog: pvacfuse valid_alleles

.. _pvacvector_valid_algorithms:

List Valid Algorithms
---------------------

.. program-output:: pvacvector valid_algorithms -h

.. .. argparse::
    :module: lib.valid_algorithms
    :func: define_parser
    :prog: pvacvector valid_algorithms

List Allele-Specific Cutoffs
----------------------------

.. program-output:: pvacvector allele_specific_cutoffs -h
