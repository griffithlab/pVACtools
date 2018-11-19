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
   contact
   mailing_list

New in version |release|
------------------------

This is a hotfix release. It fixes the following issue(s):

- When using the MHCnuggets prediction algorithm for MHC class II alleles
  (``MHCnuggetsII``) not all epitope sequences were predicted for inframe
  insertions. This issues has now been fixed.
- For MHCflurry cases with peptide sequences that were shorter than the
  desired epitope length were not handled correctly which resulted in an
  error. This issues has been resolved in this release.

New in version 1.1.2
------------------------

This is a hotfix release. It fixes the following issue(s):

- In version 1.1.0 we added a ``--pass-only`` flag to pVACseq that would
  result in only variants with ``FILTER`` of ``PASS`` or ``.`` getting processed.
  However, this option was not getting passed along to the pVACseq process
  correctly, resulting in this option not taking effect. This hotfix release
  fixes this issue and the ``--pass-only`` flag should now work as expected.

New in version 1.1.1
--------------------

This is a hotfix release. It fixes the following issue(s):

- In version 1.1 we updated VAFs to be fractions, rather than percentages. A
  bug in this code change resulted in an error when using custom VAF cutoff
  values instead of the default. This has now been fixed.

New in version |version|
------------------------

This version adds a host of new features to pVACtools:

- pVACseq is now able to parse VAF, depth, and expression information directly
  from the VCF. This makes the ``--additional-input-file-list`` option
  obsolete. The ``--additional-input-file-list`` option is now deprecated and will be removed in an
  upcoming release. For more information on how to annotate your VCF with
  readcount and expression information, see the :ref:`prerequisites_label` page.
- pVACseq is now able to handle proximal germline and somatic variants. In
  order to incorporate those into the epitope predictions, you will need to
  provide a phased variants VCF to your pVACseq run using the
  ``--phased-proximal-variants-vcf`` option. For more information on how to
  create this file, see the :ref:`prerequisites_label` page.
- We added support to pVACseq for filtering on transcript support levels. This requires
  the input VCF to be annotated with the TSL field by VEP. Be default, any
  transcripts with a TSL above 1 will be filtered out.
- The binding filter of pVACseq and pVACfuse can now be run with flexible, allele-specific
  binding-thresholds. This feature can be enabled using the
  ``--allele-specific-binding-thresholds`` flag. The thresholds used are taken
  from `the IEDB recommendations
  <https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions>`_.
- pVACseq now supports a ``--pass-only`` flag that will result in any VCF
  entries with a ``FILTER`` to be skipped. Using this flag, only VCF entries
  with a ``FILTER`` of ``PASS`` or ``.`` will be processed.
- We added support for the `MHCflurry <http://www.biorxiv.org/content/early/2017/08/09/174243>`_ and
  `MHCnuggets <http://karchinlab.org/apps/appMHCnuggets.html>`_ prediction algorithms. These
  can be used by listing ``MHCflurry``, ``MHCnuggetsI`` (for MHC Class I alleles),
  and/or ``MHCnuggetsII`` (for MHC Class II alleles) as the prediction
  algorithms in your run commands.
- The default ``--tdna-vaf`` and ``--trna-vaf`` cutoff values have been
  updated from 0.4 to 0.25. This is the minimum VAF threshold that an epitope
  candidate must meet in order to pass the coverage filter.
- We now offer a graphical user interface, :ref:`pvacviz`, to run pVACseq as an alernative
  to using the command line. pVACviz, can also be used to plot and filter your pVACseq
  results.

To stay up-to-date on the latest pVACtools releases please join our :ref:`mailing_list`.

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
