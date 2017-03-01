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

This version improves the sorting of the final report file. The file will now be
sorted by "Gene Name" and "Mutation" and within these categories by the MT score
- either "Median MT Score" or "Best MT score" depending on the
top-score-metric used.

This release implements a 60 second wait between each request to the IEDB
RESTful API in order to decrease the load on their servers. We recommend the
usage of the :ref:`standalone IEDB tools <iedb_install>` for long-running processes with many
variants, prediction algorithm, epitope lengths, or alleles.

There also have been various bugfixes in this release.

Citation
--------

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632 <http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
