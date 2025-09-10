pVACtools
=========

pVACtools is a cancer immunotherapy tools suite consisting of the following
tools:

**pVACseq**
   A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a VCF file.

**pVACbind**
   A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a FASTA file.

**pVACfuse**
   A tool for detecting neoantigens resulting from gene fusions.

**pVACsplice**
   A tool for detecting neoantigens resulting from splice site variants.

**pVACvector**
   A tool designed to aid specifically in the construction of DNA-based
   cancer vaccines.

**pVACview**
   An application based on R Shiny that assists
   users in reviewing, exploring and prioritizing neoantigens from the results of
   pVACtools processes for personalized cancer vaccine design.

.. image:: images/pVACtools_main-figure_v8.png
    :align: center
    :alt: pVACtools immunotherapy workflow

Contents
--------

.. toctree::
   :maxdepth: 2

   pvacseq
   pvacbind
   pvacfuse
   pvacsplice
   pvacvector
   pvacview

.. toctree::
   :maxdepth: 1

   install
   courses
   tools
   frequently_asked_questions
   releases
   license
   citation
   funding
   contribute
   contact
   mailing_list

New in Version |release|
------------------------

This is a bugfix release. It fixes the following problem(s):

* Fix a couple of issues with the new `--top-score-metric2` by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1291

  * When adding the ``--top-score-metric2`` option, the logic for determining
    the included candidates during aggregate report creation was amended to
    compare either the IC50 or percentile to the aggregate inclusion binding
    threshold. This logic should not have been changed and instead only the
    IC50 should be compared to the aggregate inclusion binding threshold, no
    matter which ``--top-score-metric2`` was selected. This specific change
    has been reverted
  * In order to achieve deterministic results when using the percentile
    ``--top-score-metric2`` option, a peptides.sort() call was used. This
    release replaces this with a better way of finding the best peptide by
    using the IC50 as a secondary sort criteria.

New in Version 5
----------------

This is a major version release. Please note that pVACtools 5.0 is not guaranteed to be
backwards-compatible and certain changes could break old workflows.

New Tools
_________

This release adds a new tool, pVACsplice, for prediction neoantigens from
splice sites. Please see the :ref:`full tool documentation <pvacsplice>` for more information.
by @mrichters in https://github.com/griffithlab/pVACtools/pull/911

New Features
____________

- This release refactors the pVACvector graph building algorithm in order to increase the probability
  for finding a solution and reduce the number of iterations needed before a solution is found. Please
  see the PR describtion for the full details. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1163
- Add a new ``--aggregate-inclusion-count-limit`` parameter to set the maximum number of epitopes
  to include in the metrics.json detailed data for variants that have a large number of candidate
  neoantigens (e.g., frameshifts). by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1147
- Add a new ``--biotypes`` parameter which defines a list of transcript consequence biotypes that the
  predictions from pVACseq and pVACsplice should be limited to. by @mrichters in https://github.com/griffithlab/pVACtools/pull/911
- Add support for additional alleles that weren't previously supported, includings ones for dog,
  mouse, and MHC class II. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1148

Bugfixes
________

- This relase fixes a bug with the ``--agregate-inclusion-binding-threshold`` which would not be used if
  the ``--allele-specific-binding-thresholds`` flag was set. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1147
- The pVACview percentile plots have been updated to include percentiles from elution and immunogenicity
  algorithms. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1149
- This release fixes a bug where the incorrect neoantigen fasta entry may be used for the reference proteome
  search if there were multiple variants or alt alleles located at the same genomic position. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1153
- Add additional trailing amino acids for frameshift insertions when creating fasta in order to capture a
  matched wildtype entry in large repetitive regions. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1155

Past release notes can be found on our :ref:`releases` page.

To stay up-to-date on the latest pVACtools releases please join our :ref:`mailing_list`.

Citations
---------

Jasreet Hundal , Susanna Kiwala , Joshua McMichael, Chris Miller, Huiming Xia,
Alex Wollam, Conner Liu, Sidi Zhao, Yang-Yang Feng, Aaron Graubert, Amber Wollam,
Jonas Neichin, Megan Neveau, Jason Walker, William Gillanders,
Elaine Mardis, Obi Griffith, Malachi Griffith. pVACtools: A Computational Toolkit to
Identify and Visualize Cancer Neoantigens. Cancer Immunology Research.
2020 Mar;8(3):409-420. doi: 10.1158/2326-6066.CIR-19-0401.
PMID: `31907209 <https://www.ncbi.nlm.nih.gov/pubmed/31907209>`_.

Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, William C. Chapman,
Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith.
`Accounting for proximal variants improves neoantigen prediction <https://www.nature.com/articles/s41588-018-0283-9>`_.
Nature Genetics. 2018, DOI: 10.1038/s41588-018-0283-9. PMID: `30510237 <https://www.ncbi.nlm.nih.gov/pubmed/30510237>`_.

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi
L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVACseq: A genome-guided
in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016,
8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632
<http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

Huiming Xia, My H. Hoang, Evelyn Schmidt, Susanna Kiwala, Joshua McMichael, Zachary L. Skidmore, Bryan Fisk, Jonathan J. Song, Jasreet Hundal, Thomas Mooney, Jason R. Walker, S. Peter Goedegebuure, Christopher A. Miller, William E. Gillanders, Obi L. Griffith,  Malachi Griffith. `pVACview: an interactive visualization tool for efficient neoantigen prioritization and selection <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01384-7>`_. Genome Medicine. 2024, 16:132, DOI: 10.1186/s13073-024-01384-7. PMID: `39538339 <http://www.ncbi.nlm.nih.gov/pubmed/39538339>`_. 

Source code
-----------
The pVACtools source code is available in `GitHub <https://github.com/griffithlab/pVACtools>`_.

License
-------
This project is licensed under `BSD 3-Clause Clear License <https://spdx.org/licenses/BSD-3-Clause-Clear.html>`_.
