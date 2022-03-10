.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

.. _pvacfuse_tools:

Optional Downstream Analysis Tools
==================================

Generate Protein Fasta
----------------------

.. program-output:: pvacfuse generate_protein_fasta -h

This tool will extract protein sequences surrounding fusion variant in an by parsing Integrate-Neo or AGFusion
output. One use case for this tool is to help select long peptides that contain short neoepitope 
candidates. For example, if pvacfuse was run to predict nonamers (9-mers) that are good binders and
the user wishes to select long peptide (e.g. 24-mer) sequences that contain the nonamer for synthesis
or encoding in a DNA vector. The fusion position will be centered in the protein sequence returned (if possible).
If the fusion causes a frameshift, the full downstream protein sequence will be returned unless the user specifies otherwise 
as described above.

Generate Aggregated Report
--------------------------

.. program-output:: pvacfuse generate_aggregated_report -h

This tool produces an aggregated version of the all_epitopes TSV. It finds the best-scoring
(lowest binding affinity) epitope for each variant, and outputs additional binding affinity for that epitope.
It also gives information about the total number of well-scoring epitopes for each variant,
as well as the HLA alleles that those epitopes are well-binding to.
For a full overview of the output, see the pVACfuse :ref:`output file documentation <pvacfuse_aggregated>`.

Calculate Reference Proteome Similarity
---------------------------------------

.. program-output:: pvacfuse calculate_reference_proteome_similarity -h

This tool will Blast peptides against the relative reference proteome and return the results in an output
TSV & reference_match file pair, given a pVACfuse run's fasta and filtered/all_epitopes TSV.  Typically,
this can be done as part of the pVACfuse run pipeline for the filtered output TSV if specified.  This tool,
however, provides a standalone way to run this on pVACfuse's generated filtered/all_epitopes TSV files.  For
instance, this may be desired if pvacfuse was originally run without this specified and one wished to
perform this additional step after the fact for the filtered TSV—or perhaps instead the results of this were
desired for the all_epitopes TSV which does not have this step performed.
For a closer look at the generated reference_match file,
see the pVACfuse :ref:`output file documentation <pvacfuse_reference_matches>`.

NetChop Predict Cleavage Sites
------------------------------

.. program-output:: pvacfuse net_chop -h

This tool uses NetChop to predict cleavage sites for neoepitopes from a pVACfuse run's filtered/all_epitopes
TSV.  In its output, it adds to the TSV 3 columns: Best Cleavage Position, Best Cleavage Score, and a
Cleavage Sites list.  Typically this step is done in the pVACfuse run pipeline for the filtered output TSV
when specified.  This tool provides a way to manually run this on pVACfuse's generated filtered/all_epitopes
TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACfuse in
the :ref:`output file documentation <pvacfuse_all_ep_and_filtered>`.

NetMHCStab Predict Stability
----------------------------

.. program-output:: pvacfuse netmhc_stab -h

This tool uses NetMHCstabpan to add stability predictions for neoepitopes from a pVACfuse run's
filtered/all_epitopes TSV.  In its output, it adds to the TSV 4 columns: Predicted Stability, Half Life,
Stability Rank, and NetMHCStab Allele.  Typically this step is done in the pVACfuse run pipeline for the
filtered output TSV when specified.  This tool provides a way to manually run this on pVACfuse's generated
filtered/all_epitopes TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACfuse in
the :ref:`output file documentation <pvacfuse_all_ep_and_filtered>`.
