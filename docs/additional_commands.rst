Additional Commands
===================

To make using pVAC-Seq easier several convenience methods are included in the package.

.. _example_data:

Download Example Data
---------------------

.. argparse::
    :module: lib.download_example_data
    :func: define_parser
    :prog: pvacseq download_example_data

.. _install_vep_plugin_label:

Install VEP Plugin
------------------

.. argparse::
    :module: lib.install_vep_plugin
    :func: define_parser
    :prog: pvacseq install_vep_plugin

List Valid Alleles
------------------

.. argparse::
    :module: lib.valid_alleles
    :func: define_parser
    :prog: pvacseq valid_alleles

Documentation For Configuration Files
-------------------------------------

.. argparse::
    :module: lib.config_files
    :func: define_parser
    :prog: pvacseq config_files
