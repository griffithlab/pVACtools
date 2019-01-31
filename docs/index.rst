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

New in version |version|
------------------------

This version adds a few features and updates:

- pVACvector now accepts a list of spacers to use when testing junction
  epitopes. These can be specified using the ``--spacers`` parameter with a
  comma-separated list of spacer peptides. Including the string ``None`` will
  also test each junction without spacers. The default is
  ``None,HH,HHC,HHH,HHHD,HHHC,AAY,HHHH,HHAA,HHL,AAL``
- The ``--expn-val`` cutoff parameter has been updated to be a float instead
  of an integer. This allows the user to provide a decimal
  cuttoff for the filtering on gene and transcript expression values.
  Previously, only whole numbers were accepted.
- Decimal numbers in the pVACseq reports are now rounded to three decimal
  places. Previously, they were not rounded.

In addition, this version also fixes a few bugs:

- The ``--normal-vaf`` cutoff value was incorrectly defaulting to 0.2 instead
  of 0.02. This resulted in the coverage filter not being as stringent as it
  should've been.
- There were a number of bugs in pVACapi and pVACviz that would prevent a user
  from submitting jobs using the interface in certain conditions. These have been resolved.
- pVACseq would previsouly not support SVs in the input VCF where the alt had
  a value of ``<DEL>``. These kinds of variants are now supported.

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
