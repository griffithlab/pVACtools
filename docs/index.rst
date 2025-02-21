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

This is a minor feature release. It adds the following features:

- Add a new parameter ``--percentile-threshold-strategy`` that controls
  filtering and tiering behavior when a percentile threshold is set.
  If this parameter is set to ``conservative`` a candidate has to pass both
  the binding affinity and percentile score thresholds. If it is set to
  ``exploratory`` it only has to pass either the binding affinity or the
  percentile score threshold. by @ldhtnp in https://github.com/griffithlab/pVACtools/pull/1185
- Add a new parameter ``--netmhciipan-version`` that controls which version
  of ``NetMHCIIpan`` and ``NetMHCIIpanEL`` are being run. The default remains
  version 4.1. by @ldhtnp in https://github.com/griffithlab/pVACtools/pull/1181
- The meaning of the Mutation Position column (in the all_epitopes.tsv and filtered.tsv pVACseq reports)
  and the Pos column (in the aggregated pVACseq report) has been updated to reflect the position(s) in
  the mutant epitope that are different from the matched wildtype epitope.

  - This is different from the previous behaviors, particular for indels, where previously
    this column was trying to reflect where the mutation occurred. However, this has proven
    difficult to programmatically determine correctly for cases with proximal variants,
    in repetitive regions, or where the mutant amino acid(s) are similar the wildtype amino
    acid(s).
  - For inframe indels, this change might now result in some positions getting marked as different
    from the matched wildtype amino acid even though they aren't technically
    part of the variant because wildtype amino acids are shifted in relation to the mutant in this
    type of variants. These shifted amino acids, although not mutated, are now different
    between mutant and matched wildtype epitope and marked as such. We do believe that this is the better approach
    because it allows us to evaluate the absolute differences between the matched wildtype and mutant epitopes.
  - These columns are now always ``NA`` if the wildtype epitope is ``NA``
  - For making the anchor position evaluation, epitopes with more than two Mutation Position/Pos entries are auto-passed.

It also fixes the following problem(s):

- Fix a bug in the Mutation Position/Pos column where proximal variants and
  certain inframe indels might result in the incorrect position(s) being
  identified. As part of this bugfix, the meaning of these columns has been
  updated (see above.)
- Fix Mutation Position column by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1206
- Fix pVACvector bug that would result in not all junctional epitopes getting tested on clipping. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1200
- Fix a pVACvector bug that would result in an error if the ``-k`` flag was not set in a run. by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1186
- Fix warning related to new lines in the mouse anchor data files in pVACview. by @ldhtnp in https://github.com/griffithlab/pVACtools/pull/1184
- Fix the column tooltips in pVACview which were previously shifted by a couple of columns. by @ldhtnp in https://github.com/griffithlab/pVACtools/pull/1189
- Fix description of the ``ref_fasta`` positional parameter in pvacsplice. by @mhoang22 in https://github.com/griffithlab/pVACtools/pull/1201
- Fix an issue in pVACsplice which would throw an error with non-human transcripts. by @mhoang22 in https://github.com/griffithlab/pVACtools/pull/1198

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
