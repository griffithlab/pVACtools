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
Additionally, variants can be limited to specific transcript biotypes
using the ``--biotypes`` parameters, which is set to only include ``protein_coding``
transcripts by default.

The output can be further limited to only certain variants by providing
a pVACseq report file to the ``--input-tsv`` argument. Only the peptide sequences for the epitopes in the TSV
will be used when creating the FASTA. If this argument is an aggregated TSV
file, use the ``--aggregate-report-evaluation`` parameter to only include
peptide sequences for epitopes matching the chosen Evaluation(s). This is
useful when creating a peptide fasta for vaccine ordering after using pVACview
to select vaccine candidates and exporting the results to a new TSV.

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

Add ML Predictions
------------------

.. program-output:: pvacseq add_ml_predictions -h

This tool adds machine learning (ML)-based neoantigen prioritization predictions to existing pVACseq output files. 
It uses a trained random forest model to predict whether neoantigen candidates should be evaluated as "Accept", 
"Reject", or "Pending" based on a comprehensive set of features derived from binding affinity predictions, 
expression data, and variant characteristics.

This tool requires that you have already generated both MHC Class I and Class II aggregated reports using 
the ``generate_aggregated_report`` command or by running the pVACseq pipeline (``pvacseq run``). It takes as input 
the Class I aggregated TSV, Class I all epitopes TSV, and Class II aggregated TSV files from a pVACseq run. 
The tool merges these files, performs data cleaning and imputation, and applies the ML model to generate evaluation predictions for each variant.

The output file is named ``<sample_name>_predict_pvacview.tsv`` and contains all columns from the original 
Class I aggregated file with two additional columns:  


.. list-table::

 * - ``Evaluation``
   - The ML-predicted evaluation status: "Accept", "Reject", or "Pending", based on the prediction probability score.
 * - ``ML Prediction (score)``
   - A formatted string combining the model-predicted evaluation with the prediction probability score (e.g., 
     "Accept (0.72)"). It shows "NA" for variants where the model could not make a prediction, which may be due to a candidate 
     not being present in either the Class I or Class II aggregated reports.

The ``--threshold_accept`` parameter controls the probability threshold for Accept predictions (default: 0.55). 
Variants with prediction probabilities >= this threshold are evaluated as "Accept". The ``--threshold_reject`` parameter 
controls the probability threshold for Reject predictions (default: 0.30). Variants with prediction probabilities <= 
this threshold are evaluated as "Reject". Everything in between is set to "Pending" for manual review. 
The ``--artifacts_path`` parameter allows you to specify a custom directory 
containing ML model artifacts, though by default the tool uses the model artifacts included with the pvactools package.

Calculate Reference Proteome Similarity
---------------------------------------

.. program-output:: pvacseq calculate_reference_proteome_similarity -h

This tool will find matches of the epitope candidates in the reference proteome and return the results in an output
TSV & reference_match file pair. It requires the input of a pVACseq run's fasta file in order to look up the larger
peptide sequence the epitope was derived from. Any substring of that peptide
sequence that matches against the reference proteome and is at least as long as the specified match length, will be
considered a hit. This tool also requires the user to provide a filtered.tsv,
all_epitopes.tsv or aggregated.tsv pVACseq report file as an input and any
candidates in this input file will be searched for.

This tool may be either run with BLASTp using either the ``refseq_select_prot`` or ``refseq_protein`` database.
By default this option uses the BLAST API but users may :ref:`independently install BLASTp <blast>`. Alternatively, users
may provide a reference proteome fasta file and this tool will string match on
the entries of this fasta file directly. This approach is recommended, because
it is significantly faster than BLASTp. Reference proteome fasta files may be
downloaded from Ensembl. For example, the latest reference proteome fasta for human
can be downloaded from `this
link <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz>`_.

For more details on the  generated reference_match file,
see the pVACseq :ref:`output file documentation <reference_matches>`.

NetChop Predict Cleavage Sites
------------------------------

.. program-output:: pvacseq net_chop -h

This tool uses NetChop to predict cleavage sites for neoepitopes from a pVACseq run's filtered/all_epitopes
TSV.  In its output, it adds to the TSV 3 columns: Best Cleavage Position, Best Cleavage Score, and a
Cleavage Sites list.  Typically this step is done in the pVACseq run pipeline for the filtered output TSV
when specified.  This tool provides a way to manually run this on pVACseq's generated filtered/all_epitopes
TSV files so that you can add this information when not present if desired.

You can view more information about these columns for pVACseq in
the :ref:`output file documentation <all_ep_and_filtered>`.

NetMHCStab Predict Stability
----------------------------

.. program-output:: pvacseq netmhc_stab -h

This tool uses NetMHCstabpan to add stability predictions for neoepitopes from a pVACseq run's
filtered/all_epitopes TSV.  In its output, it adds to the TSV 4 columns: Predicted Stability, Half Life,
Stability Rank, and NetMHCStab Allele.  Typically this step is done in the pVACseq run pipeline for the
filtered output TSV when specified.  This tool provides a way to manually run this on pVACseq's generated
filtered/all_epitopes TSV files so that you can add this information when not present if desired.

You can view more information about these columns for pVACseq in
the :ref:`output file documentation <all_ep_and_filtered>`.

Identify Problematic Amino Acids
--------------------------------

.. program-output:: pvacseq identify_problematic_amino_acids -h

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

This tool may be used with any filtered.tsv or all_epitopes.tsv pVACseq report
file.

Mark Genes of Interest
----------------------

.. program-output:: pvacseq mark_genes_of_interest -h

Update Tiers
------------

.. program-output:: pvacseq update_tiers -h

Create Peptide Ordering Form
----------------------------

.. program-output:: pvacseq create_peptide_ordering_form -h

This tool generates a comprehensive peptide ordering package from pVACseq results in a single step. 
It streamlines the preparation of long peptides for synthesis by combining protein sequence extraction, 
manufacturability assessment, peptide annotation, and visualization into one workflow. The output includes 
peptide FASTA files, manufacturability reports, and color-coded Excel summaries that highlight binding strength, 
sequence properties, and variant context.

This command replaces the need to run the ``generate_protein_fasta``, ``generate_reviews_files``,
and ``color_peptides51mer`` scripts separately. The output includes the following files:

.. list-table::

 * - ``<output_file>_<sample_name>.fa``
   - Contains the generated peptides in FASTA format for peptide synthesis.
 * - ``<output_file>_<sample_name>.manufacturability.tsv``
   - Manufacturability assessments for the peptides, including metrics such as cysteine
     content, hydrophobicity, and sequence complexity.
 * - ``<output_file>_<sample_name>.Colored_Peptides.xlsx``
   - A color-coded Excel file summarizing peptides, annotations, manufacturability metrics,
     and peptide positions, ready for ordering.
 * - ``<output_file>_<sample_name>.Annotated.Neoantigen_Candidates.xlsx``
   - A spreadsheet intended for downstream manual review of selected variants, including
     visualization in tools such as IGV.

Several options are available for tailoring the output. The ``flanking_sequence_length`` determines the number of flanking amino acids around the mutation of interest when creating the peptide sequence for the ordering spreadsheet. The ``--biotypes`` option can be used to pre-filter transcripts when generating peptide sequences
from the input VCF, limiting the analysis to specific transcript biotypes (default: ``protein_coding``) and the ``--pass-only`` flag can be used to narrow down the input VCF to PASS variants only. The two latter options should match the options selected for the original pVACseq run so that variants will match correctly between the peptide sequences created by this tool and the variants in the classI_tsv aggregated report.

Additionally, the ``--aggregate-report-evaluation`` parameter can be used
to restrict the output reports to candidates with specific evaluation states in the classI_tsv (e.g. ``Accept``, ``Reject``,
``Pending``, or ``Review``; multiple values may be provided as a comma-separated list). 

For custom peptide prioritization thresholds, the IC50 and percentile cutoffs for class I and II
can be adjusted using the appropriate flags.
