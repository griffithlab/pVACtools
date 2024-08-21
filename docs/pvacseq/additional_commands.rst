.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Additional Commands
===================

To make using pVACseq easier, several convenience methods are included in the package.

.. _example_data:

Download Example Data
---------------------

.. program-output:: pvacseq download_example_data -h

.. .. argparse::
    :module: lib.download_example_data
    :func: define_parser
    :prog: pvacseq download_example_data

.. _install_vep_plugin_label:

Install VEP Plugin
------------------

.. program-output:: pvacseq install_vep_plugin -h

.. .. argparse::
    :module: lib.install_vep_plugin
    :func: define_parser
    :prog: pvacseq install_vep_plugin

.. _valid_alleles:

List Valid Alleles
------------------

.. program-output:: pvacseq valid_alleles -h

.. .. argparse::
    :module: lib.valid_alleles
    :func: define_parser
    :prog: pvacseq valid_alleles

.. _valid_algorithms:

List Valid Algorithms
---------------------

.. program-output:: pvacseq valid_algorithms -h

.. .. argparse::
    :module: lib.valid_algorithms
    :func: define_parser
    :prog: pvacseq valid_algorithms

List Allele-Specific Cutoffs
----------------------------

.. program-output:: pvacseq allele_specific_cutoffs -h
