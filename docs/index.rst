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

**pVACviz**
   A browser-based user interface that assists
   users in launching, managing, reviewing, and visualizing the results of
   pVACtools processes.

.. image:: images/pVACtools_main-figure_v3d.png
    :align: center
    :alt: pVACtools immunotherapy workflow


.. toctree::
   :maxdepth: 2

   pvacseq
   pvacfuse
   pvacvector
   pvacviz

.. toctree::
   :maxdepth: 1

   install
   frequently_asked_questions
   releases
   citation
   contact
   mailing_list

New in version |release|
------------------------

This version is a hotfix release. It fixes the following issues:

- Tensorflow is incompatible with multiprocessing when the parent process
  imports tensorflow or a tensorflow-dependent module. For this reason
  MHCflurry and MHCnuggets were removed from parallelization. In this
  release we moved to calling MHCflurry and MHCnuggets on the command line,
  which allowed us to remove our direct imports of these modules and allows us
  to parallelize the calls to these two prediction algorithms. All prediction
  algorithms supported by pVACtools can now be used in multiprocessing mode.
- Some users were reporting ``Illegal instruction (core dumped)`` errors
  because their hardware was incompatible with the version of tensorflow we
  were using. Pinning the tensorflow version to 1.5.0 with this release should
  solve this problem.
- When running in multiprocessing mode while using the IEDB API, users would
  experience a higher probability of failed requests to the API. The IEDB API
  would throw a 403 error when rejecting requests due to too
  many simultaneous requests. pVACtools would previously not retry on this type of
  error. This release now adds retries on this error code. We also improved
  the random wait time calculation between requests so that the likelihood of
  multiple retries hitting at the same time has now been reduced.
- When encountering a truncated input VCF, the VCF parser used by pVACtools
  would throw an error that was not indicative of the real error source.
  pVACseq now catches these errors and emmits a more descriptive error message
  when encountering a truncated VCF.
- One option when annotating a VCF with VEP is the ``-total-length`` flag. When
  using this flag, the total length would be written to the
  ``Protein_position`` field. pVACseq previously did not support a VCF with a
  ``Protein_position`` field in this format. This release adds support for it.
- When creating the combined MHC class I and MHC class II all_epitopes file,
  we were previously not correctly determining all necessary headers which
  would lead to incorrect output of the individual prediction algorithm score
  columns. This release fixes this issue.

Past release notes can be found on our :ref:`releases` page.

To stay up-to-date on the latest pVACtools releases please join our :ref:`mailing_list`.

Citations
---------

Jasreet Hundal, Susanna Kiwala, Joshua McMichael, Christopher A Miller,
Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng,
Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker,
William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith.
`pVACtools: a computational toolkit to select and visualize cancer
neoantigens <https://doi.org/10.1101/501817>`_.
bioRxiv 501817; doi: https://doi.org/10.1101/501817

Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, William C. Chapman, Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. `Accounting for proximal variants improves neoantigen prediction <https://www.nature.com/articles/s41588-018-0283-9>`_. Nature Genetics. 2018, DOI: 10.1038/s41588-018-0283-9. PMID: `30510237 <https://www.ncbi.nlm.nih.gov/pubmed/30510237>`_.

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi
L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVACseq: A genome-guided
in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016,
8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632
<http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
