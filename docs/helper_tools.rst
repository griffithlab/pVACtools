.. _helper_tools:

Helper Tools
============

This section describes a set of utility tools designed to assist users in working more efficiently with the core functionality of the software.
These helper tools provide streamlined workflows, support data formatting and validation, and simplify common post-processing steps. Our comparison
tool and others included here can help users compare results, generate reports, and prepare data for further analysis.

.. toctree::
   :glob:

   helper_tools/comparison_tool.rst

.. _download_wdls:

Download WDL tool wrappers
--------------------------

Workflow Description Language (WDL) tool wrappers for pVACseq and pVACfuse
can be downloaded using the ``pvactools download_wdls`` command.

.. program-output:: pvactools download_wdls -h

.. _download_cwls:

Download CWL tool wrappers
--------------------------

Common Workflow Language (CWL) tool wrappers for pVACseq, pVACfuse, and
pVACvector can be downloaded using the ``pvactools download_cwls`` command.

.. program-output:: pvactools download_cwls -h

.. _valid_alleles:

List Valid Alleles
------------------

.. program-output:: pvactools valid_alleles -h

.. _valid_algorithms:

List Valid Algorithms
---------------------

.. program-output:: pvactools valid_algorithms -h

List Valid NetMHCIIpan/NetMHCIIpanEL Versions
---------------------------------------------

.. program-output:: pvactools valid_netmhciipan_versions -h

List Allele-Specific Cutoffs
----------------------------

.. program-output:: pvactools allele_specific_cutoffs -h
