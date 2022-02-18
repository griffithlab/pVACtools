.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

.. _pvacbind_tools:

Optional Downstream Analysis Tools
==================================

Generate Aggregated Report
--------------------------

.. program-output:: pvacbind generate_aggregated_report -h

This tool produces an aggregated version of the all_epitopes TSV. It finds the best-scoring
(lowest binding affinity) epitope for each variant, and outputs additional binding affinity for that epitope.
It also gives information about the total number of well-scoring epitopes for each variant,
as well as the HLA alleles that those epitopes are well-binding to.
For a full overview of the output, see the pVACbind :ref:`output file documentation <pvacbind_aggregated>`.

Calculate Reference Proteome Similarity
---------------------------------------

.. program-output:: pvacbind calculate_reference_proteome_similarity -h

This tool will Blast peptides against the relative reference proteome and return the results in an output
TSV & reference_match file pair, given a pVACbind run's fasta and filtered/all_epitopes TSV.  Typically,
this can be done as part of the pVACbind run pipeline for the filtered output TSV if specified.  This tool,
however, provides a standalone way to run this on pVACbind's generated filtered/all_epitopes TSV files.  For
instance, this may be desired if pvacbind was originally run without this specified and one wished to
perform this additional step after the fact for the filtered TSV—or perhaps instead the results of this were
desired for the all_epitopes TSV which does not have this step performed.
For a closer look at the generated reference_match file,
see the pVACbind :ref:`output file documentation <pvacbind_reference_matches>`.

NetChop Predict Cleavage Sites
------------------------------

.. program-output:: pvacbind net_chop -h

This tool uses NetChop to predict cleavage sites for neoepitopes from a pVACbind run's filtered/all_epitopes
TSV.  In its output, it adds to the TSV 3 columns: Best Cleavage Position, Best Cleavage Score, and a
Cleavage Sites list.  Typically this step is done in the pVACbind run pipeline for the filtered output TSV
when specified.  This tool provides a way to manually run this on pVACbind's generated filtered/all_epitopes
TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACbind in
the :ref:`output file documentation <pvacbind_all_ep_and_filtered>`.

NetMHCStab Predict Stability
----------------------------

.. program-output:: pvacbind netmhc_stab -h

This tool uses NetMHCstabpan to add stability predictions for neoepitopes from a pVACbind run's
filtered/all_epitopes TSV.  In its output, it adds to the TSV 4 columns: Predicted Stability, Half Life,
Stability Rank, and NetMHCStab Allele.  Typically this step is done in the pVACbind run pipeline for the
filtered output TSV when specified.  This tool provides a way to manually run this on pVACbind's generated
filtered/all_epitopes TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACbind in
the :ref:`output file documentation <pvacbind_all_ep_and_filtered>`.
