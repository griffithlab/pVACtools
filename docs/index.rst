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

**pVACvector**
   A tool designed to aid specifically in the construction of DNA-based
   cancer vaccines.

**pVACview**
   An application based on R Shiny that assists
   users in reviewing, exploring and prioritizing neoantigens from the results of
   pVACtools processes for personalized cancer vaccine design.

.. image:: images/pVACtools_main-figure_v7.png
    :align: center
    :alt: pVACtools immunotherapy workflow

Contents
--------

.. toctree::
   :maxdepth: 2

   pvacseq
   pvacbind
   pvacfuse
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

New in Release |release|
------------------------

This is a minor feature release. It adds the following features:

- Add a new helper command ``pvacseq|pvacfuse|pvacbind|pvacvector valid_algorithms``
  by @ldhtnp in https://github.com/griffithlab/pVACtools/pull/1108
- When running the ``pvacseq generate_protein_fasta`` command with the ``--phased-proximal-variants-vcf``
  argument, output the intermediate ``proximal_variants.tsv`` file by @evelyn-schmidt
  in https://github.com/griffithlab/pVACtools/pull/1091
- In pVACview, clear the comment text input box after saving the comment by @ldhtnp
  in https://github.com/griffithlab/pVACtools/pull/1113
- Add support for mouse allele anchor positions by @ldhtnp in
  https://github.com/griffithlab/pVACtools/pull/1110
- Skip variants where VEP didn't predict an amino acid change by @susannasiebert
  in https://github.com/griffithlab/pVACtools/pull/1121
- Update the ordering of the fasta file output of the ``pvacseq|pvacfuse generate_protein_fasta``
  command when running with the ``--input-tsv`` argument so that the order of the fasta sequences
  is consistent with the order of the neoantigen candidates in the input TSV by @mhoang22 in
  https://github.com/griffithlab/pVACtools/pull/1002
- Updat the ``pvacfuse generate_protein_fasta`` command to allow aggregated TSVs as an input TSV
  by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1134
- Update pVACview to display the anchor positions currently applied to the data by @susannasiebert
  in https://github.com/griffithlab/pVACtools/pull/1114

This release also fixes the following bug(s):

- Handle invalid pVACfuse characters by trimming the sequence instead of skipping it. The previous
  implementation would lead to missing sequences in certain downstream steps, resulting in errors.
  by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1130
- Add new pVACview R files to the list of files getting copied into the pVACseq output folder.
  These files were previsouly not copied in the the results folder, leading to error when running
  the ``pvacview run`` commands on a pVACseq output directory. by @susannasiebert in
  https://github.com/griffithlab/pVACtools/pull/1126
- Remove single DP and DQ alpha and beta chain alleles from the list of supported alleles in MHCnuggetsII.
  This is because those alleles need to be defined as a pair of alpha- and beta-chains in order to be
  meaningful. Also remove DRA alleles from the same list since the DR locus is defined only by the beta
  chain because functional variation in mature DRA gene products is absent. by @susannasiebert in
  https://github.com/griffithlab/pVACtools/pull/1133
- Fix errors in the rounding of the min and max values of the sliders in the custom pVACview module by
  @evelyn-schmidt in https://github.com/griffithlab/pVACtools/pull/1116
- Remove unused code in the Frameshift.pm VEP plugin that causes errors with certain types of variants
  by @susannasiebert in https://github.com/griffithlab/pVACtools/pull/1122

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

Source code
-----------
The pVACtools source code is available in `GitHub <https://github.com/griffithlab/pVACtools>`_.

License
-------
This project is licensed under `BSD 3-Clause Clear License <https://spdx.org/licenses/BSD-3-Clause-Clear.html>`_.
