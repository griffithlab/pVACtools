pVACtools
=========

pVACtools is a cancer immunotherapy tools suite consisting of the following
tools:

**pVACseq**
   A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a list of tumor mutations.

**pVACfuse**
   A tool for detecting neoantigens resulting from gene fusions.

**pVACvector**
   A tool designed to aid specifically in the construction of DNA-based
   cancer vaccines.

.. image:: images/pVACtools_main-figure_v2e.png
    :align: center
    :alt: pVACtools immunotherapy workflow


.. toctree::
   :maxdepth: 2

   pvacseq
   pvacfuse
   pvacvector

.. toctree::
   :maxdepth: 1

   install
   frequently_asked_questions
   contact

New in version |release|
------------------------

This is a hotfix release. It fixes the following issues:

- The epitope length used for generating the peptide fasta when running with
  multiple epitope lengths was incorrect. This would potentially result in including
  fasta sequences that were shorter than the largest epitope length which
  would cause an error during calls to IEDB.
- pVACseq would fail with a nondescript error message if the input VCF was not
  annotated with VEP before running. A more descriptive error message has been
  added.
- IEDB changed the format of class II IEDB alleles which would cause an error
  when running with those alleles. pVACtools will now handle transposing the
  affected alleles into the new format.
- The standalone binding filters had a few bugs that would result in syntax
  errors during runtime.
- The indexes created for each fusion entry with pVACfuse had the potential to
  not be unique which would result in parsing errors downstream.
- pVACseq had the potential to use the incorrect VEP allele for positions with
  multiple alternate alleles which would result in the incorrect CSQ entry
  getting used for some of those alternate alleles.
- pVACseq would throw an error if the chosen peptide sequence length exceeds
  the wildtype protein sequence length of a transcript.

Coming soon
-----------

**pVACclient**
   A browser-based user interface that assists
   users in launching, managing, reviewing, and visualizing the results of
   pVACtools processes.

**pVACapi**
    The pVACapi will provide a HTTP REST interface to the pVACtools
    suite.

Citation
--------

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi
L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVACseq: A genome-guided
in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016,
8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632
<http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
