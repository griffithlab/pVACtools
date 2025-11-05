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

Identify Problematic Amino Acids
--------------------------------

.. program-output:: pvacbind identify_problematic_amino_acids -h

This tool is used to identify positions in an epitope with an amino acid that
is problematic for downstream processing, e.g. vaccine manufacturing. Since
this can differ from case to case, this tool requires the user to specify which
amino acid(s) to consider problematic. This can be specified in one of three
formats:

.. list-table::

 * - ``amino_acid(s)``
   - One or more one-letter amino acid codes. Any occurrence of this amino acid string,
     regardless of the position in the epitope, is problematic. When specifying more than
     one amino acid, they will need to occur together in the specified order.
 * - ``amino_acid:position``
   - A one letter amino acid code, followed by a colon separator, followed by a positive
     integer position (one-based). The occurrence of this amino acid at the position
     specified is problematic., E.g. G:2 would check for a Glycine at the second position
     of the epitope. The N-terminus is defined as position 1.
 * - ``amino_acid:-position``
   - A one letter amino acid code, followed by a colon separator, followed by a negative
     integer position. The occurrence of this amino acid at the specified position from
     the end of the epitope is problematic. E.g., G:-3 would check for a Glycine at the
     third position from the end of the epitope. The C-terminus is defined as position -1.

You may specify any number of these problematic amino acid(s), in any
combination, by providing them as a comma-separated list.

This tool may be used with any filtered.tsv or all_epitopes.tsv pVACbind report
file.

Update Tiers
------------

.. program-output:: pvacseq update_tiers -h
