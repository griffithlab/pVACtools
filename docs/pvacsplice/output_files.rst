.. image:: ../images/pVACsplice_logo_trans-bg_v4b.png
    :align: right
    :alt: pVACsplice logo
    :width: 175px

.. _pvacsplice_output_files:

Output Files
============

The pVACsplice pipeline will write a few files to the main output directory.
These files contain information about the processed splice sites that aren't
specific to either class of prediction algorithm:

- ``<sample_name>.transcripts.fa`` and matching ``.fai`` index file: A fasta file of wildtype and splice site peptide sequences for splice sites predicted by RegtTools and supported by pVACsplice.
- ``<sample_name>_combined.tsv``: A TSV file combining information for each splice site from the RegTools TSV and the input VCF.

pVACsplice writes its results in separate folders depending on
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
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with (optionally) cleavage site, and stability
       predictions added.
   * - ``<sample_name>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format. Not
       generated when running with elution algorithms only.
   * - ``<sample_name>.all_epitopes.aggregated.tsv.reference_matches`` (optional)
     - A file outlining details of reference proteome matches

Additionally, each folder will contain subfolders, one for each selected
epitope length, that contains intermediate files that are specific to each
epitope length.

Filters applied to the filtered.tsv file
----------------------------------------

The filtered.tsv file is the all_epitopes file with the following filters
applied (in order):

- Binding Filter
- Coverage Filter
- Transcript Support Level Filter
- Top Score Filter

Please see the :ref:`Standalone Filter Commands<pvacsplice_filter_commands>`
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

Prediction Algorithms Supporting Elution Scores
_______________________________________________

- MHCflurryEL
- NetMHCpanEL
- NetMHCIIpanEL
- BigMHC_EL

Prediction Algorithms Supporting Immunogenicity Scores
______________________________________________________

- BigMHC_IM
- DeepImmuno

Please note that when running pVACseq with only elution or immunogenicity algorithms, no
aggregate report and pVACview files are created.

.. _pvacsplice_all_ep_and_filtered:

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
   * - ``Junction``
     - Junction ID in regtools output
   * - ``Junction Start``
     - The start position of this junction in the zero-based, half-open coordinate system
   * - ``Junction Stop``
     - The stop position of this junction in the zero-based, half-open coordinate system
   * - ``Junction Score``
     - The number of reads supporting the junction. (integer)
   * - ``Junction Anchor``
     - Field that specifies the donor, acceptor configuration. See `Notes <https://regtools.readthedocs.io/en/latest/commands/junctions-annotate/#notes>`_ (D/A/DA/NDA/N) 
   * - ``Transcript``
     - The Ensembl ID of the affected transcript
   * - ``Transcript Support Level``
     - The `transcript support level (TSL) <https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl>`_
       of the affected transcript. ``Not Supported`` if the VCF entry doesn't contain TSL information.
   * - ``Biotype``
     - The biotype of the affected transcript
   * - ``Ensembl Gene ID``
     - The Ensembl ID of the affected gene
   * - ``Variant Type``
     - The type of variant. ``missense`` for missense mutations, ``inframe_ins`` for
       inframe insertions, ``inframe_del`` for inframe deletions, and ``FS`` for frameshift variants
   * - ``Amino Acid Change``
     - The amnio acid change of this mutation
   * - ``Gene Name``
     - The Ensembl gene name of the affected gene
   * - ``HGVSc``
     - The HGVS coding sequence variant name
   * - ``HGVSp``
     - The HGVS protein sequence variant name
   * - ``WT Protein Length``
     - Length of fully-translated wildtype protein
   * - ``ALT Protein Length``
     - Length of fully-translated alternate protein
   * - ``Frameshift Event``
     - Is the variant a frameshift event? (yes/no)
   * - ``Protein Position``
     - Starting position of Epitope (Position of the first amino acid of selected epitope in the fully-translated protein)
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Peptide Length``
     - The peptide length of the epitope
   * - ``Epitope Seq``
     - The mutant epitope sequence
   * - ``Median IC50 Score``
     - Median ic50 binding affinity of the mutant epitope across all prediction algorithms used
   * - ``Best IC50 Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Best IC50 Score Method``
     - Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
   * - ``Median Percentile``
     - Median binding affinity percentile rank of the mutant epitope across all prediction algorithms (those that provide percentile output)
   * - ``Best Percentile``
     - Lowest percentile rank of this epitope's ic50 binding affinity of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile Method``
     - Prediction algorithm with the lowest binding affinity percentile rank for this epitope
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
   * - ``Index``
     - A unique idenitifer for this variant-transcript combination
   * - ``Fasta Key``
     - the number identifier for corresponding altered peptide isoform in pvac output fasta
   * - ``Individual Prediction Algorithm  IC50 Scores and Percentiles`` (multiple)
     - ic50 binding affintity and percentile ranks for the ``Epitope Seq`` for the individual prediction algorithms used
   * - ``MHCflurryEL WT and MT Processing Score and Presentation Score and Percentile`` (optional)
     - MHCflurry elution processing score and presentation score and percentiles
       for the ``Epitope Seq`` if the run included
       MHCflurryEL as one of the prediction algorithms
   * - ``Problematic Positions`` (optional)
     - A list of positions in the ``MT Epitope Seq`` that match the
       problematic amino acids defined by the ``--problematic-amino-acids``
       parameter
   * - ``Gene of Interest`` (T/F)
     - Is the ``Gene Name`` found in the genes of interest list?
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

.. _pvacsplice_aggregated:

all_epitopes.aggregated.tsv Report Columns
--------------------------------------------

The ``all_epitopes.aggregated.tsv`` file is an aggregated version of the all_epitopes TSV.
It shows the :ref:`best-scoring epitope <pvacsplice_best_peptide>`
for each variant, and outputs additional binding affinity, expression, and
coverage information for that epitope. It also gives information about the
total number of well-scoring epitopes for each variant, the number of
transcripts covered by those epitopes, as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.

Only epitopes meeting the ``--aggregate-inclusion-binding-threshold`` are included in this report (default: 5000).
If the number of unique epitopes for a variant meeting this threshold exceeds the
``--aggregate-inclusion-count-limit``, only the n best-binding epitopes up to this
limit are included (default: 15). If the Best Peptide does not meet the aggregate inclusion criteria, it will be still be
counted in the ``Num Included Peptides``.

Whether the median or the lowest binding affinity metrics are used for determining the
included eptiopes, selecting the best-scoring epitope, and which values are output in the ``IC50 MT``
and ``%ile MT`` columns is controlled by the ``--top-score-metric`` parameter.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description

   * - ``ID``
     - A unique identifier for the junction (Gene name . transcript. Junction ID . variant chr . variant start - variant stop . junction type)
   * - ``HLA Alleles`` (multiple)
     - For each HLA allele in the run, the number of this variant's epitopes that bound well
       to the HLA allele (with median/lowest mutant binding affinity < binding_threshold)
   * - ``Gene``
     - The Ensembl gene name of the affected gene
   * - ``Transcript``
     - The Ensembl ID of the affected transcript
   * - ``Junction Name``
     - junction ID from regtools output
   * - ``AA Change``
     - The amino acid change for the mutation
   * - ``Best Peptide``
     - The best-binding mutant epitope sequence (see Best Peptide Criteria
       below for more details on how this is determined)
   * - ``TSL``
     - The Transcript Support Level of the Transcript
   * - ``Allele``
     - The Allele that the Best Peptide is binding to
   * - ``Pos``
     - The one-based position of the start of the mutation within the epitope sequence. ``0`` if the
       start of the mutation is before the epitope (as can occur downstream of frameshift mutations)
   * - ``Prob Pos``
     - A list of positions in the Best Peptide that are problematic.
       ``None`` if the ``--problematic-pos`` parameter was not set during
       the pVACseq run
   * - ``Num Included Peptides``
     - The number of included peptides according to the
       ``--aggregate-inclusion-binding-threshold`` and
       ``--aggregate-inclusion-count-limit``
   * - ``Num Passing Peptides``
     - The number of unique well-binding peptides for this mutation.
   * - ``IC50 MT``
     - Median or lowest ic50 binding affinity of the best-binding mutant epitope across all prediction algorithms used
   * - ``%ile MT``
     - Median or lowest binding affinity percentile rank of the best-binding mutant epitope across all prediction algorithms used (those that provide percentile output)
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
   * - ``Ref Match`` (T/F) (optional)
     - Was there a match of the mutated peptide sequence to the reference proteome?
   * - ``Evaluation``
     - Column to store the evaluation of each variant when evaluating the run in pVACview. Either ``Accept``, ``Reject``, or ``Review``.

.. _pvacsplice_best_peptide:

Best Peptide Criteria
_____________________

To determine the Best Peptide, all peptides meeting the
``--aggregate-inclusion-threshold`` and ``--aggregate-inclusion-count-limit``
(see above) are evaluated as follows:

- Pick the entries with no Problematic Positions.
- Of the remaining entries, pick the one with the lowest median/best MT IC50
  score.

.. _pvacsplice_aggregate_report_tiers_label:

The pVACsplice Aggregate Report Tiers
_____________________________________

Tiering Parameters
******************

To tier the Best Peptide, several cutoffs can be adjusted using arguments provided to the pVACsplice run:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``--binding-threshold``
     - The threshold used for filtering epitopes on the IC50 MT binding affinity.
     - 500
   * - ``--allele-specific-binding-thresholds``
     - Instead of the hard cutoff set by the ``--binding-threshold``, use
       allele-specific binding thresholds. For alleles where no
       allele-specific binding threshold is available, use the
       ``--binding-threshold`` as a fallback. To print a list of alleles that have
       specific binding thresholds and the value of those thresholds, run ``pvacseq allele_specific_cutoffs``.
     - False
   * - ``--percentile-threshold``
     - When set, use this threshold to filter epitopes on the %ile MT score in addition to having to meet the binding threshold.
     - None
   * - ``--percentile-threshold-strategy``
     - Specify the candidate inclusion strategy. The ``conservative`` option requires a candidate to pass BOTH the binding threshold
       and percentile threshold (if set). The ``exploratory`` option requires a candidate to pass EITHER the binding threshold or
       the percentile threshold.
     - conservative
   * - ``--tumor-purity``
     - Value between 0 and 1 indicating the fraction of tumor cells in the tumor sample. Information is used for a simple estimation of
       whether variants are subclonal or clonal based on VAF. If not provided, purity is estimated directly from the VAFs.
     - None
   * - ``--trna-vaf``
     - Tumor RNA VAF Cutoff. Used to calculate the allele expression cutoff for tiering.
     - 0.25
   * - ``--trna-cov``
     - Tumor RNA Coverage Cutoff. Used as a cutoff for tiering.
     - 10
   * - ``--expn-val``
     - Gene and Expression cutoff. Used to calculate the allele expression cutoff for tiering.
     - 1.0
   * - ``--maximum-transcript-support-level``
     - The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL).
       Transcript support level needs to be <= this cutoff to be included in most tiers.
     - 1

Tiers
*****

Given the thresholds provided above, the Best Peptide is evaluated and binned into a tier as follows:

.. list-table::
   :header-rows: 1

   * - Tier
     - Citeria
   * - ``Pass``
     - Best Peptide passes the binding, reference match, expression, transcript, clonal, and problematic position criteria
   * - ``PoorBinder``
     - Best Peptide fails the binding criteria but passes the reference match, expression, transcript, clonal, and problematic position criteria
   * - ``RefMatch``
     - Best Peptide fails the reference match criteria but passes the binding, expression, transcript, clonal, and problematic position criteria
   * - ``PoorTranscript``
     - Best Peptide fails the transcript criteria but passes the binding, reference match, expression, clonal, and problematic position criteria
   * - ``LowExpr``
     - Best Peptide meets the low expression criteria and passes the binding, reference match, transcript, clonal, and problematic position criteria
   * - ``Subclonal``
     - Best Peptide fails the clonal criteria but passes the binding, reference match, expression, transcript, and problematic position criteria
   * - ``ProbPos``
     - Best Peptide fails the problematic position criteria but passes the binding, reference match, expression, transcript, and clonal criteria
   * - ``Poor``
     - Best Peptide doesn't fit in any of the above tiers, usually if it fails
       two or more criteria
   * - ``NoExpr``
     - Best Peptide is not expressed (RNA Expr == 0 or RNA VAF == 0)


Criteria Details
****************

.. list-table::
   :header-rows: 1

   * - Criteria
     - Description
     - Evaluation Logic
   * - Binding Criteria
     - Pass if Best Peptide is strong binder
     - ``IC50 MT < binding_threshold`` and ``%ile MT < percentile_threshold``
       (if ``--percentile-threshold`` parameter is set and 'conservative' ``--percentile-threshold-strategy`` is used) or
       ``IC50 MT < binding_threshold`` or ``%ile MT < percentile_threshold``
       (if 'exploratory' ``--percentile-threshold-strategy`` is used)
   * - Expression Criteria
     - Pass if Best Transcript is expressed
     - ``Allele Expr > trna_vaf * expn_val``
   * - Reference Match Criteria
     - Pass if there are no reference protome matches
     - ``Ref Match == True``
   * - Transcript Criteria
     - Pass if Best Transcript matches any of the user-specified ``--transcript-prioritization-strategy`` criteria
     - ``TSL <= maximum_transcript_support_level`` (if
       ``--transcript-prioritization-strategy`` includes ``tsl``)
       ``MANE Select == True`` (if ``--transcript-prioritization-strategy
       includes ``mane_select``)
       ``Canonical == True`` (if ``--transcript-prioritization-strategy``
       incluces ``canonical``)
   * - Low Expression Criteria
     - Peptide has low expression or no expression but RNA VAF and coverage
     - ``(0 < Allele Expr < trna_vaf * expn_val) OR (RNA Expr == 0 AND RNA
       Depth > trna_cov AND RNA VAF > trna_vaf)``
   * - Clonal Criteria
     - Best Peptide is likely in the founding clone of the tumor
     - ``DNA VAF > tumor_purity / 4``
   * - Problematic Position Criteria
     - Best Peptide contains a problematic amino acid as defined by the
       ``--problematic-amino-acids`` parameters
     - ``Prob Pos == None``

.. _pvacsplice_reference_matches:

aggregated.tsv.reference_matches Report Columns
-----------------------------------------------

This file is only generated when the ``--run-reference-proteome-similarity``
option is chosen.

.. flat-table::
   :header-rows: 1

   * - Column Name
     - Description (BLAST)
     - Description (reference fasta)
   * - ``Chromosome``
     - :cspan:`2` The chromosome of this variant
   * - ``Start``
     - :cspan:`2` The start position of this variant in the zero-based, half-open coordinate system
   * - ``Stop``
     - :cspan:`2` The stop position of this variant in the zero-based, half-open coordinate system
   * - ``Reference``
     - :cspan:`2` The reference allele
   * - ``Variant``
     - :cspan:`2` The alt allele
   * - ``Transcript``
     - :cspan:`2` The Ensembl ID of the affected transcript
   * - ``MT Epitope Seq``
     - :cspan:`2` The mutant peptide sequence for the epitope candidate
   * - ``Peptide``
     - The peptide sequence submitted to BLAST
     - The peptide sequence to search for in the reference proteome
   * - ``Hit ID``
     - The BLAST alignment hit ID (reference proteome sequence ID)
     - The FASTA header ID of the entry where the match was made
   * - ``Hit Definition``
     - The BLAST alignment hit definition (reference proteome sequence name)
     - The FASTA header description of the entry where the match was made
   * - ``Match Window``
     - :cspan:`2` The substring of the ``Peptide`` that was found in the ``Match
       Sequence``
   * - ``Match Sequence``
     - The BLAST match sequence
     - The FASTA sequence of the entry where the match was made
   * - ``Match Start``
     - :cspan:`2` The match start position of the ``Match Window`` in the ``Match Sequence``
   * - ``Match Stop``
     - :cspan:`2` The match stop position of the ``Match Window`` in the ``Match Sequence``

