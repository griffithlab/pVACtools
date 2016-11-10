.. pVAC-Seq documentation master file, created by
   sphinx-quickstart on Wed Aug 31 08:31:01 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pVAC-Seq
====================================

pVAC-Seq is a cancer immunotherapy pipeline for the identification of **p**\ ersonalized **V**\ ariant **A**\ ntigens by **C**\ ancer **Seq**\ uencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. It is being used in studies of checkpoint therapy response and to identify targets for cancer vaccines and adoptive T cell therapies. For more general information, see the `manuscript published in Genome Medicine <http://www.genomemedicine.com/content/8/1/11>`_.

.. toctree::
   :maxdepth: 3

   features
   install
   prerequisites
   run
   filter_commands
   additional_commands
   optional_downstream_analysis_tools
   contact

New in version |version|
------------------------

pVAC-Seq now supports local installs of IEDB MHC `class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_ binding prediction tools. This feature can be used by passing the directory that contains the local installations to the ``--iedb-install-directory`` parameter.

This version adds a new column ``Mutation Position`` to the report output. This column denotes the 1-based start position of the mutation in the ``MT Epitope Seq``. If the value is ``0`` the mutation start position is before the first position in the epitope.

pVAC-Seq now allows the user to specify the number of retries after a request to the IEDB RESTful interface fails. The number of retries can be set by using the ``--iedb-retries`` parameter. Previously this number was hard-coded to 3. More retries might be necessary in order to get a successful response for complex queries (e.g., large number of variants, long frameshift downstream sequences, choice of compute-intensive prediction algorithms). This parameter should be used in conjunction with ``--fasta-size`` and ``--downstream-sequence-length`` for the highest likelihood of success of finishing a pVAC-Seq run.

This release fixes an error that was introduced in the previous version which would occur when the user would try to rerun a process in the same output directory.

This version also fixes a bug with how to handle variants that are no-call or homozygous-reference. These variants will now be skipped.

Citation
--------

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632 <http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
