.. image:: images/pVACview_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACview logo

.. _pvacview:

pVACview
========

pVACview is a R shiny based tool designed to aid specifically in the prioritization and selection of neoantigen candidates for personalized cancer vaccines.
It takes as inputs a pVACseq output aggregate report file (tsv format) and a corresponding pVACseq output metrics file (json).
pVACview allows the user to launch an R shiny application to load and visualize the given neoantigen candidates with detailed information
including that of the genomic variant, transcripts covering the variant, and good-binding peptides predicted from the respective transcripts.
It also incorporates anchor prediction data for a range of class I HLA alleles and peptides ranging from 8 to 11-mers.
By taking all levels of information into account for the neoantigen candidates, clinicians will be able to make more informed
decisions when deciding final peptide candidates for personalized cancer vaccines.


.. toctree::
   :maxdepth: 2
   :glob:

   pvacview/prerequisites
   pvacview/pvacseq_module
   pvacview/neofox_module
   pvacview/custom_module
   pvacview/troubleshooting
