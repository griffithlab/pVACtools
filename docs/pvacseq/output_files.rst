.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

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
       proccessable variant-transcript combinations.
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with (optionally) cleavage site, stability
       predictions, and reference proteome similarity metrics added.
   * - ``<sample_name>.filtered.tsv.reference_matches`` (optional)
     - A file outlining details of reference proteome matches
   * - ``<sample_name>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format.

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
     - The one-based position of the start of the mutation within the epitope sequence. ``0`` if the
       start of the mutation is before the epitope
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
     - Prediction algorithm with the lowest ic50 binding affinity percentile rank for this epitope
   * - ``Best MT Percentile``
     - Lowest percentile rank of this epitope's ic50 binding affinity of all prediction algorithms used
   * - ``Corresponding WT Percentile``
     - ic50 binding affinity percentile rank of the wildtype epitope. ``NA`` if there is no ``WT Epitope Seq``.
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
     - Median ic50 binding affinity percentile rank of the mutant epitope across all prediction algorithms used
   * - ``Median WT Percentile``
     - Median ic50 binding affinity percentile rank of the wildtype epitope across all prediction algorithms used.
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
     -  The match start position in the matched reference proteome sequence
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
   * - ``HLA Alleles`` (multiple) (T/F)
     - For each HLA allele in the run, did the mutation result in an epitope that binded well
       to the HLA allele (median mutant binding affinity < 1000).
   * - ``Gene``
     - The Ensembl gene name of the affected gene
   * - ``AA_change``
     - The amino acid change for the mutation
   * - ``Num_Transcript``
     - The Number of transcripts for this mutation that resulted in at least
       one well-binding peptide (median mutant binding affinity < 1000).
   * - ``Peptide``
     - The best-binding mutant epitope sequence (lowest median mutant binding
       affinity)
   * - ``Pos``
     - The one-based position of the start of the mutation within the epitope sequence. ``0`` if the
       start of the mutation is before the epitope
   * - ``Num_Peptides``
     - The number of unique well-binding peptides for this mutation.
   * - ``ic50_MT``
     - Median ic50 binding affinity of the best-binding mutant epitope across all prediction algorithms used
   * - ``ic50_WT``
     - Median ic50 binding affinity of the corresponding wildtype epitope across all prediction algorithms used.
   * - ``percentile_MT``
     - Median ic50 binding affinity percentile rank of the best-binding mutant epitope across all prediction algorithms used
   * - ``percentile_WT``
     - Median ic50 binding affinity percentile rank of the corresponding wildtype epitope across all prediction algorithms used.
   * - ``RNA_expr``
     - Gene expression value for the annotated gene containing the variant.
   * - ``RNA_VAF``
     - Tumor RNA variant allele frequency (VAF) at this position.
   * - ``RNA_Depth``
     - Tumor RNA depth at this position.
   * - ``DNA_VAF``
     - Tumor DNA variant allele frequency (VAF) at this position.
   * - ``tier``
     - A tier corresponding to the suitability of variants for use in vaccines.

The pVACseq Aggregate Report Tiers
__________________________________

To bin a variant in a tier, the best binding epitope is evaluated as follows:

.. list-table::
   :header-rows: 1

   * - Tier
     - Citeria
   * - ``Pass``
     - Median MT Score < 500 and (Tumor RNA VAF * Gene Expression) > 3 and
       Tumor DNA VAF > (clonal VAF / 2) and not an anchor residue
   * - ``Relaxed``
     - Median MT Score < 1000 and (Tumor RNA VAF * Gene Expression) > 1 and
       Tumor DNA VAF > (clonal VAF / 2) and not an anchor residue
   * - ``Anchor``
     - Median MT Score < 1000 and (Tumor RNA VAF * Gene Expression) > 1 and
       Tumor DNA VAF > (clonal VAF / 2) and an anchor residue
   * - ``Subclonal``
     - Median MT Score < 1000 and (Tumor RNA VAF * Gene Expression) > 1 and
       Tumor DNA VAF < (clonal VAF / 2) and not an anchor residue
   * - ``LowExpr``
     - Median MT Score < 1000 and low expression and
       Tumor DNA VAF > (clonal VAF / 2) and not an anchor residue
   * - ``NoExpr``
     - Gene Expression == 0 and Tumor DNA VAF == 0 and expression is not low
   * - ``Poor``
     - Variant did not match any of the above criteria

.. list-table::
   :header-rows: 0

   * - anchor residue
     - Anchor residue positions are definied as the first, second, second-to-last, and
       last position in the epitope. The mutation in an epitope is an anchor residue if the mutation position is at one
       of those positions and the median wildtype binding affinity is < 1000.
   * - low expression
     - (Tumor RNA VAF * Gene Expression) > 0 or (Gene Expression == 0 and Tumor
       RNA Depth < 50 and Tumor RNA VAF > .10)
