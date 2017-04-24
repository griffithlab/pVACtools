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

This release adds handling for DNPs and MNPs missense mutations.

This version adds a new option ``--additonal-report-columns`` to the ``pvacseq
run`` command which can be use
to append additional columns of data to the report. Right now the only value
supported for this option is ``sample_name`` which appends a column with the
sample name to the final report.

We updated the logic that determines whether or not a corresponding wildtype
epitope for a mutant epitope is included in the report. Previously, we would only
include the corresponding wildtype epitope if the number of **consecutive**
matching amino acids between mutant and wildtype epitope was larger then half
of the total number of amino acids in the epitope. We now use the **total** number of
matching amino acids between the mutant epitope and the corespondig wildtype epitope
across the whole length of the epitope to make that determination. The total
number of matching amino acids needs to be larger than half of the length of
the epitope. Otherwise the corresponding wildtype epitope is reported as "NA".

With this release any execution of ``pvacseq run`` will create a log file of the
inputs used. This log file is then used when executing another run
with the same output directory. This ensures that you can only write to the same
output directory if identical parameters are used.

Citation
--------

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632 <http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
