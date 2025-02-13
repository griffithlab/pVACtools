.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

Additional Commands
===================

To make using pVACbind easier, several convenience methods are included in the package.

.. _pvacbind_example_data:

Download Example Data
---------------------

.. program-output:: pvacbind download_example_data -h

List Valid Alleles
------------------

.. program-output:: pvacbind valid_alleles -h

.. _pvacbind_valid_algorithms:

List Valid Algorithms
---------------------

.. program-output:: pvacbind valid_algorithms -h

.. .. argparse::
    :module: lib.valid_algorithms
    :func: define_parser
    :prog: pvacbind valid_algorithms

List Valid NetMHCIIpan/NetMHCIIpanEL Versions
---------------------------------------------

.. program-output:: pvacbind valid_netmhciipan_versions -h

List Allele-Specific Cutoffs
----------------------------

.. program-output:: pvacbind allele_specific_cutoffs -h
