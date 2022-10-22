.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _pvacseq_output_files:

Output Files
============

The pVACseq pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms
- ``combined``: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoepitope predictions from both

Each folder will contain the same list of output files (listed in the order
created):

.. list-table::
   :header-rows: 1

   * - File Name
     - Description
   * - ``<sample_name>.tsv``
     - An intermediate file with variant, transcript, coverage, vaf, and expression
       information parsed from the input files.
   * - ``<sample_name>.tsv_<chunks>`` (multiple)
     - The above file but split into smaller chunks for easier processing with IEDB.
   * - ``<sample_name>.fasta``
     - A fasta file with mutant and wildtype peptide subsequences for all
       processable variant-transcript combinations.
   * - ``<sample_name>.net_chop.fa``
     - A fasta file with mutant and wildtype peptide subsequences specific for use in running the net_chop tool.
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``. Only
       epitopes resulting from supported variants (missense, inframe indels, and frameshifts)
       are included. If the ``--pass-only`` flag is
       set, variants that have a FILTER set in the VCF are excluded.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with (optionally) cleavage site, stability
       predictions, and reference proteome similarity metrics added.
   * - ``<sample_name>.filtered.tsv.reference_matches`` (optional)
     - A file outlining details of reference proteome matches
   * - ``<sample_name>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format.
   * - ``ui.R``, ``app.R``, ``server.R``, ``styling.R``, ``anchor_and_helper_functions.R``
     - pVACview R Shiny application files
   * - ``www`` (directory)
     - Directory containing image files for pVACview

Filters applied to the filtered.tsv file
----------------------------------------

The filtered.tsv file is the all_epitopes file with the following filters
applied (in order):

- Binding Filter
- Coverage Filter
- Transcript Support Level Filter
- Top Score Filter

Please see the :ref:`Standalone Filter Commands<filter_commands>`
documentation for more information on each individual filter. The standalone
filter commands may be useful to reproduce the filtering or to chose different
filtering thresholds.

Prediction Algorithms Supporting Percentile Information
_______________________________________________________

pVACseq outputs binding affinity percentile rank information when provided by
a chosen prediction algorithm. The following prediction algorithms calculate a
percentile rank:

- MHCflurry
- NetMHC
- NetMHCcons
- NetMHCpan
- NetMHCIIpan
- NNalign
- PickPocket
- SMM
- SMMPMBEC
- SMMalign

The following prediction algorithms do not provide a percentile rank:

- MHCnuggets

.. _all_ep_and_filtered:

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Chromosome``
     - The chromosome of this variant
   * - ``Start``
     - The start position of this variant in the zero-based, half-open coordinate system
   * - ``Stop``
     - The stop position of this variant in the zero-based, half-open coordinate system
   * - ``Reference``
     - The reference allele
   * - ``Variant``
     - The alt allele
   * - ``Transcript``
     - The Ensembl ID of the affected transcript
   * - ``Transcript Support Level``
     - The `transcript support level (TSL) <https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl>`_
       of the affected transcript. ``NA`` if the VCF entry doesn't contain TSL information.
   * - ``Ensembl Gene ID``
     - The Ensembl ID of the affected gene
   * - ``Variant Type``
     - The type of variant. ``missense`` for missense mutations, ``inframe_ins`` for
       inframe insertions, ``inframe_del`` for inframe deletions, and ``FS`` for frameshift variants
   * - ``Mutation``
     - The amnio acid change of this mutation
   * - ``Protein Position``
     - The protein position of the mutation
   * - ``Gene Name``
     - The Ensembl gene name of the affected gene
   * - ``HGVSc``
     - The HGVS coding sequence variant name
   * - ``HGVSp``
     - The HGVS protein sequence variant name
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Peptide Length``
     - The peptide length of the epitope
   * - ``Sub-peptide Position``
     - The one-based position of the epitope within the protein sequence used to make the prediction
   * - ``Mutation Position``
     - The one-based positional range (inclusive) of the mutation within the epitope sequence.  If the mutation is a deletion, the amino acids flanking the deletion are recorded. ``0`` represents that some or all of the mutation is before the epitope, ``length+1`` represents some or all of the mutation is after the epitope, otherwise it indexes specific amino acid(s) within the epitope.  Note that in the case of ambiguous amino acid changes, this reflects the change that is left-aligned, starting from the first changed amino acid; this may differ from the ``Mutation`` column.
   * - ``MT Epitope Seq``
     - The mutant epitope sequence
   * - ``WT Epitope Seq``
     - The wildtype (reference) epitope sequence at the same position in the full protein sequence. ``NA`` if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated
   * - ``Best MT Score Method``
     - Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
   * - ``Best MT Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Corresponding WT Score``
     - ic50 binding affinity of the wildtype epitope. ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Corresponding Fold Change``
     - ``Corresponding WT Score`` / ``Best MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Best MT Percentile Method``
     - Prediction algorithm with the lowest binding affinity percentile rank for this epitope
   * - ``Best MT Percentile``
     - Lowest percentile rank of this epitope's ic50 binding affinity of all prediction algorithms used (those that provide percentile output)
   * - ``Corresponding WT Percentile``
     - binding affinity percentile rank of the wildtype epitope. ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Tumor DNA Depth``
     - Tumor DNA depth at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
   * - ``Tumor DNA VAF``
     - Tumor DNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain
       tumor DNA readcount annotation.
   * - ``Tumor RNA Depth``
     - Tumor RNA depth at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
   * - ``Tumor RNA VAF``
     - Tumor RNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain
       tumor RNA readcount annotation.
   * - ``Normal Depth``
     - Normal DNA depth at this position. ``NA`` if VCF entry does not contain normal DNA readcount annotation.
   * - ``Normal VAF``
     - Normal DNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain
       normal DNA readcount annotation.
   * - ``Gene Expression``
     - Gene expression value for the annotated gene containing the variant. ``NA`` if VCF entry does not contain
       gene expression annotation.
   * - ``Transcript Expression``
     - Transcript expression value for the annotated transcript containing the variant. ``NA`` if VCF entry does
       not contain transcript expression annotation.
   * - ``Median MT Score``
     - Median ic50 binding affinity of the mutant epitope across all prediction algorithms used
   * - ``Median WT Score``
     - Median ic50 binding affinity of the wildtype epitope across all prediction algorithms used.
       ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Median Fold Change``
     - ``Median WT Score`` / ``Median MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Median MT Percentile``
     - Median binding affinity percentile rank of the mutant epitope across all prediction algorithms (those that provide percentile output)
   * - ``Median WT Percentile``
     - Median binding affinity percentile rank of the wildtype epitope across all prediction algorithms used (those that provide percentile output)
       ``NA`` if there is no ``WT Epitope Seq``.
   * - ``Individual Prediction Algorithm WT and MT Scores and Percentiles`` (multiple)
     - ic50 binding affintity and percentile ranks for the ``MT Epitope Seq`` and ``WT Eptiope Seq`` for the individual prediction algorithms used
   * - ``Index``
     - A unique idenitifer for this variant-transcript combination
   * - ``cterm_7mer_gravy_score``
     - Mean hydropathy of last 7 residues on the C-terminus of the peptide
   * - ``max_7mer_gravy_score``
     - Max GRAVY score of any kmer in the amino acid sequence. Used to determine if there are any extremely
       hydrophobic regions within a longer amino acid sequence.
   * - ``difficult_n_terminal_residue`` (T/F)
     - Is N-terminal amino acid a Glutamine, Glutamic acid, or Cysteine?
   * - ``c_terminal_cysteine`` (T/F)
     - Is the C-terminal amino acid a Cysteine?
   * - ``c_terminal_proline`` (T/F)
     - Is the C-terminal amino acid a Proline?
   * - ``cysteine_count``
     - Number of Cysteines in the amino acid sequence. Problematic because they can form disulfide bonds across
       distant parts of the peptide
   * - ``n_terminal_asparagine`` (T/F)
     - Is the N-terminal amino acid a Asparagine?
   * - ``asparagine_proline_bond_count``
     - Number of Asparagine-Proline bonds. Problematic because they can spontaneously cleave the peptide
   * - ``Best Cleavage Position`` (optional)
     - Position of the highest predicted cleavage score
   * - ``Best Cleavage Score`` (optional)
     - Highest predicted cleavage score
   * - ``Cleavage Sites`` (optional)
     - List of all cleavage positions and their cleavage score
   * - ``Predicted Stability`` (optional)
     - Stability of the pMHC-I complex
   * - ``Half Life`` (optional)
     - Half-life of the pMHC-I complex
   * - ``Stability Rank`` (optional)
     - The % rank stability of the pMHC-I complex
   * - ``NetMHCstab allele`` (optional)
     - Nearest neighbor to the ``HLA Allele``. Used for NetMHCstab prediction
   * - ``Reference Match`` (T/F) (optional)
     - Was there a BLAST match of the mutated peptide sequence to the
       reference proteome?

.. image:: ../images/output_file_columns.png
    :alt: pVACseq ouput file columns illustration

.. _reference_matches:

filtered.tsv.reference_matches Report Columns
---------------------------------------------

This file is only generated when the ``--run-reference-proteome-similarity``
option is chosen.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Chromosome``
     - The chromosome of this variant
   * - ``Start``
     - The start position of this variant in the zero-based, half-open coordinate system
   * - ``Stop``
     - The stop position of this variant in the zero-based, half-open coordinate system
   * - ``Reference``
     - The reference allele
   * - ``Variant``
     - The alt allele
   * - ``Transcript``
     - The Ensembl ID of the affected transcript
   * - ``Peptide``
     - The peptide sequence submitted to BLAST
   * - ``Hit ID``
     - The BLAST alignment hit ID (reference proteome sequence ID)
   * - ``Hit Definition``
     - The BLAST alignment hit definition (reference proteome sequence name)
   * - ``Query Sequence``
     - The BLAST query sequence
   * - ``Match Sequence``
     - The BLAST match sequence
   * - ``Match Start``
     - The match start position in the matched reference proteome sequence
   * - ``Match Stop``
     - The match stop position in the matched reference proteome sequence

.. _aggregated:

all_epitopes.aggregated.tsv Report Columns
--------------------------------------------

The ``all_epitopes.aggregated.tsv`` file is an aggregated version of the all_epitopes TSV.
It presents the best-scoring (lowest binding affinity)
epitope for each variant, and outputs additional binding affinity, expression, and
coverage information for that epitope. It also gives information about the
total number of well-scoring epitopes for each variant, the number of
transcripts covered by those epitopes, as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``ID``
     - A unique identifier for the variant
   * - ``HLA Alleles`` (multiple)
     - For each HLA allele in the run, the number of this variant's epitopes that bound well
       to the HLA allele (with median mutant binding affinity < 1000)
   * - ``Gene``
     - The Ensembl gene name of the affected gene
   * - ``AA Change``
     - The amino acid change for the mutation
   * - ``Num Passing Transcripts``
     - The number of transcripts for this mutation that resulted in at least
       one well-binding peptide (median mutant binding affinity < 1000).
   * - ``Best Peptide``
     - The best-binding mutant epitope sequence (lowest median mutant binding
       affinity)
   * - ``Pos``
     - The one-based position of the start of the mutation within the epitope sequence. ``0`` if the
       start of the mutation is before the epitope (as can occur downstream of frameshift mutations)
   * - ``Num Passing Peptides``
     - The number of unique well-binding peptides for this mutation.
   * - ``IC50 MT``
     - Median ic50 binding affinity of the best-binding mutant epitope across all prediction algorithms used
   * - ``IC50 WT``
     - Median ic50 binding affinity of the corresponding wildtype epitope across all prediction algorithms used.
   * - ``%ile MT``
     - Median binding affinity percentile rank of the best-binding mutant epitope across all prediction algorithms used (those that provide percentile output)
   * - ``%ile WT``
     - Median binding affinity percentile rank of the corresponding wildtype epitope across all prediction algorithms used (those that provide percentile output)
   * - ``RNA Expr``
     - Gene expression value for the annotated gene containing the variant.
   * - ``RNA VAF``
     - Tumor RNA variant allele frequency (VAF) at this position.
   * - ``Allele Expr``
     - RNA Expr * RNA VAF
   * - ``RNA Depth``
     - Tumor RNA depth at this position.
   * - ``DNA VAF``
     - Tumor DNA variant allele frequency (VAF) at this position.
   * - ``Tier``
     - A tier suggesting the suitability of variants for use in vaccines.
   * - ``Evaluation``
     - Column to store the evaluation of each variant when evaluating the run in pVACview. Either ``Accept``, ``Reject``, or ``Review``.

.. _pvacseq_aggregate_report_tiers_label:

The pVACseq Aggregate Report Tiers
__________________________________

To bin a variant in a tier, the best binding epitope is evaluated as follows:

.. list-table::
   :header-rows: 1

   * - Tier
     - Citeria
   * - ``NoExpr``
     - Mutant allele is not expressed
   * - ``LowExpr``
     - Mutant allele has low expression `(TPM * RNA_VAF < 1)`
   * - ``Subclonal``
     - Likely not in the founding clone of the tumor `(DNA_VAF > max(DNA_VAF)/2)`
   * - ``Anchor``
     - Mutation is at an anchor residue in the shown peptide, and the WT allele has good binding `(WT ic50 <1000)`
   * - ``Poor``
     - Fails two or more of the above criteria
   * - ``Relaxed``
     - Passes the above criteria, has decent MT binding `(ic50 < 1000)`
   * - ``Pass``
     - Passes the above criteria, has strong MT binding `(ic50 < 1000)` and strong expression `(TPM * RNA_VAF > 3)`
