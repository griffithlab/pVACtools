.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _optional_downstream_analysis_tools_label:

Optional Downstream Analysis Tools
==================================

Generate Protein Fasta
----------------------

.. program-output:: pvacseq generate_protein_fasta -h

.. .. argparse::
    :module: lib.generate_protein_fasta
    :func: define_parser
    :prog: pvacseq generate_protein_fasta

This tool will extract protein sequences surrounding supported protein altering variants in an
input VCF file. One use case for this tool is to help select long peptides that contain short neoepitope 
candidates. For example, if pvacseq was run to predict nonamers (9-mers) that are good binders and
the user wishes to select long peptide (e.g. 24-mer) sequences that contain the nonamer for synthesis
or encoding in a DNA vector. The protein sequence extracted will correspond to the transcript sequence 
used in the annotated VCF. The alteration in the VCF (e.g. a somatic missense SNV) will be centered in the 
protein sequence returned (if possible). If the variant is near the beginning or end of the CDS, it will
be as close to center as possible while returning the desired protein sequence length. If the variant
causes a frameshift, the full downstream protein sequence will be returned unless the user specifies otherwise 
as described above. The ``flanking_sequence_length`` positional parameter
controls how many amino acids will be included on either side of the mutation.

To incorporate proximal variants in the final sequence, use the
``--phased-proximal-variants-vcf`` argument. Please see the :ref:`phased_vcf`
section of the documentation on how to create this VCF.

The output may be limited to PASS variants only by setting the ``--pass`` only
flag and to mutant sequences by setting the ``--mutant-only`` flag.

The output can be further limited to only certain variants by providing
pVACseq report file to the ``--input-tsv`` argument.Only the peptide sequences for the epitopes in the TSV
will be used when creating the FASTA. If this argument is an aggregated TSV
file, use the ``--aggregate-report-evaluation`` parameter to only include
peptide sequences for epitopes matching the chosen Evaluation(s).

Generate Aggregated Report
--------------------------

.. program-output:: pvacseq generate_aggregated_report -h

This tool produces an aggregated version of the all_epitopes TSV. It finds the best-scoring (lowest binding affinity)
epitope for each variant, and outputs additional binding affinity, expression, and
coverage information for that epitope. It also gives information about the
total number of well-scoring epitopes for each variant, the number of
transcripts covered by those epitopes, as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.
For a full definition of these tiers, see the pVACseq :ref:`output file documentation <aggregated>`.

Calculate Reference Proteome Similarity
---------------------------------------

.. program-output:: pvacseq calculate_reference_proteome_similarity -h

This tool will Blast peptides against the relative reference proteome and return the results in an output
TSV & reference_match file pair, given a pVACseq run's fasta and filtered/all_epitopes TSV.  Typically, this
can be done as part of the pVACseq run pipeline for the filtered output TSV if specified.  This tool,
however, provides a standalone way to run this on pVACseq's generated filtered/all_epitopes TSV files.  For
instance, this may be desired if pvacseq was originally run without this specified and one wished to perform
this additional step after the fact for the filtered TSV—or perhaps instead the results of this were desired
for the all_epitopes TSV which does not have this step performed.
For a closer look at the generated reference_match file,
see the pVACseq :ref:`output file documentation <reference_matches>`.

NetChop Predict Cleavage Sites
------------------------------

.. program-output:: pvacseq net_chop -h

This tool uses NetChop to predict cleavage sites for neoepitopes from a pVACseq run's filtered/all_epitopes
TSV.  In its output, it adds to the TSV 3 columns: Best Cleavage Position, Best Cleavage Score, and a
Cleavage Sites list.  Typically this step is done in the pVACseq run pipeline for the filtered output TSV
when specified.  This tool provides a way to manually run this on pVACseq's generated filtered/all_epitopes
TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACseq in
the :ref:`output file documentation <all_ep_and_filtered>`.

NetMHCStab Predict Stability
----------------------------

.. program-output:: pvacseq netmhc_stab -h

This tool uses NetMHCstabpan to add stability predictions for neoepitopes from a pVACseq run's
filtered/all_epitopes TSV.  In its output, it adds to the TSV 4 columns: Predicted Stability, Half Life,
Stability Rank, and NetMHCStab Allele.  Typically this step is done in the pVACseq run pipeline for the
filtered output TSV when specified.  This tool provides a way to manually run this on pVACseq's generated
filtered/all_epitopes TSV files so that you can add this information when not present if desired.
You can view more about these columns for pVACseq in
the :ref:`output file documentation <all_ep_and_filtered>`.

Identify Problematic Amino Acids
--------------------------------

.. program-output:: pvacseq identify_problematic_amino_acids -h
