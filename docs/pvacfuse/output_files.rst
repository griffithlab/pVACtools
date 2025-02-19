.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

.. _pvacfuse_output_files:

Output Files
============

The pVACfuse pipeline will write its results in separate folders depending on
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
   * - ``<sample_name>.fasta``
     - A fasta file with mutant peptide subsequences for each fusion.
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with cleavage site and stability
       predictions added.
   * - ``<sample_name>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format. Not generated when running with elution algorithms only.
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
- Top Score Filter

Please see the :ref:`Standalone Filter Commands<pvacfuse_filter_commands>`
documentation for more information on each individual filter. The standalone
filter commands may be useful to reproduce the filtering or to chose different
filtering thresholds.

Prediction Algorithms Supporting Elution Scores
_______________________________________________

- MHCflurryEL (Presentation and Processing)
- NetMHCpanEL
- NetMHCIIpanEL
- BigMHC_EL

Prediction Algorithms Supporting Immunogenicity Scores
______________________________________________________

- BigMHC_IM
- DeepImmuno

Please note that when running pVACfuse with only elution or immunogenicity algorithms, no
aggregate report and pVACview files are created.

Prediction Algorithms Supporting Percentile Information
_______________________________________________________

pVACfuse outputs percentile rank information when provided by
a chosen binding affinity, elution, or immunogenicity prediction algorithm.
The following prediction algorithms calculate a
percentile rank:

- MHCflurry
- MHCflurryEL (Presentation)
- MHCnuggets
- NetMHC
- NetMHCcons
- NetMHCpan
- NetMHCpanEL
- NetMHCIIpan
- NetMHCIIpanEL
- NNalign
- PickPocket
- SMM
- SMMPMBEC
- SMMalign

.. _pvacfuse_all_ep_and_filtered:

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Chromosome``
     - The chromosomes of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Start``
     - The start positions of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Stop``
     - The stop positions of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Transcript``
     - The Ensembl IDs of the affected transcripts
   * - ``Gene Name``
     - The Ensembl gene names of the affected genes
   * - ``Variant Type``
     - The type of fusion. ``inframe_fusion`` for inframe fusions, ``frameshift_fusion`` for frameshift fusions
   * - ``Mutation``
     - A unique identifier for the fusion
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Sub-peptide Position``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Epitope Seq``
     - Epitope sequence
   * - ``Median IC50 Score``
     - Median ic50 binding affinity of the epitope of all prediction algorithms used
   * - ``Best IC50 Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Best IC50 Score Method``
     - Prediction algorithm with the lowest ic50 binding affinity for this epitope
   * - ``Median Percentile``
     - Median binding affinity percentile rank of the epitope across all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile``
     - Lowest percentile rank of this epitope's ic50 binding affinity of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile Method``
     - Prediction algorithm with the lowest binding affinity percentile rank for this epitope
   * - ``Individual Prediction Algorithm Scores and Percentiles`` (multiple)
     - ic50 binding affintity and percentile ranks for the ``Epitope Seq`` for the individual prediction algorithms used
   * - ``MHCflurryEL Processing Score and Presentation Score and Percentile`` (optional)
     - MHCflurry elution processing score and presentation score and percentiles
       for the ``Epitope Seq`` if the run included
       MHCflurryEL as one of the prediction algorithms
   * - ``Read Support``
     - The sum of spanning and encompassing reads over the fusion position.
       ``NA`` if the run was made with AGFusion data and without a
       ``--starfusion-file`` input.
   * - ``Expression``
     - The number of fusion-supporting RNA-seq fragments as FFPM (fusion fragments per million total reads). ``NA`` if the run was made
       without a ``--starfusion-file`` input.
   * - ``Problematic Positions`` (optional)
     - A list of positions in the ``Epitope Seq`` that match the
       problematic amino acids defined by the ``--problematic-amino-acids``
       parameter
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

.. _pvacfuse_aggregated:

all_epitopes.aggregated.tsv Report Columns
--------------------------------------------

The ``all_epitopes.aggregated.tsv`` file is an aggregated version of the all_epitopes TSV.
It shows the best-scoring epitope
for each variant, and outputs additional binding affinity, expression, and
coverage information for that epitope. It also gives information about the
total number of well-scoring epitopes for each variant as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.

Only epitopes meeting the ``--aggregate-inclusion-binding-threshold`` are included in this report (default: 5000).
If the number of unique epitopes for a fusion meeting this threshold exceeds the
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
     - A unique identifier for the fusion
   * - ``HLA Alleles`` (multiple)
     - For each HLA allele in the run, the number of this fusion's epitopes that bound well
       to the HLA allele (with median binding affinity < 1000)
   * - ``Gene``
     - The Ensembl gene names of the affected genes
   * - ``Best Peptide``
     - The best-binding epitope sequence (lowest ``IC50 MT`` score)
   * - ``Best Transcript``
     - The fusion transcripts coding for the Best Peptide
   * - ``Allele``
     - The Allele that the Best Peptide is binding to
   * - ``Prob Pos``
     - A list of positions in the Best Peptide that are problematic. ``None`` if the ``--problematic-pos`` parameter was not set during the pVACfuse run
   * - ``Num Included Peptides``
     - The number of included peptides according to the
       ``--aggregate-inclusion-binding-threshold`` and
       ``--aggregate-inclusion-count-limit``
   * - ``Num Passing Peptides``
     - The number of included peptides for this fisoopm that are well-binding.
   * - ``IC50 MT``
     - Median or lowest IC50 binding affinity of the best-binding epitope across all prediction algorithms used
   * - ``%ile MT``
     - Median or lowest binding affinity percentile rank of the best-binding epitope across all prediction algorithms used (those that provide percentile output)
   * - ``Expr``
     - The number of fusion-supporting RNA-seq fragments as FFPM (fusion fragments per million total reads). ``NA`` if the run was made without a ``--starfusion-file`` input.
   * - ``Read Support``
     - The sum of spanning and encompassing reads over the fusion position. ``NA`` if the run was made with AGFusion data and without a ``--starfusion-file`` input.
   * - ``Tier``
     - A tier suggesting the suitability of variants for use in vaccines.
   * - ``Ref Match`` (T/F) (optional)
     - Was there a match of the peptide sequence to the reference proteome?
   * - ``Evaluation``
     - Column to store the evaluation of each fusion. Either ``Accept``, ``Reject``, or ``Review``.

The pVACfuse Aggregate Report Tiers
___________________________________

Tiering Parameters
******************

To tier the Best Peptide, several cutoffs can be adjusted using parameters
provided to the pVACfuse run:

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
       specific binding thresholds and the value of those thresholds, run ``pvacfuse allele_specific_cutoffs``.
     - False
   * - ``--percentile-threshold``
     - When set, use this threshold to filter epitopes on the %ile MT score in addition to having to meet the binding threshold.
     - None
   * - ``--percentile-threshold-strategy``
     - Specify the candidate inclusion strategy. The ``conservative`` option requires a candidate to pass BOTH the binding threshold
       and percentile threshold (if set). The ``exploratory`` option requires a candidate to pass EITHER the binding threshold or
       the percentile threshold.
     - conservative
   * - ``--read-support``
     - The threshold used for filtering epitopes on the Read Support.
     - 5
   * - ``--expn-val``
     - The threshold used for filtering epitopes on the Expr.
     - 0.1

Tiers
*****

Given the thresholds provided above, the Best Peptide is evaluated and binned
into tiers as follows:

.. list-table::
   :header-rows: 1

   * - Tier
     - Criteria
   * - ``Pass``
     - Best Peptide passes the binding, read support, and expression criteria
   * - ``LowReadSupport``
     - Best Peptide fails the read support criteria but passes the binding and
       expression criteria
   * - ``LowExpr``
     - Best Peptide fails the expression criteria but passes the binding and
       read support criteria
   * - ``Poor``
     - Best Peptide doesn't fit any of the above tiers, usually if it fails two
       or more criteria or if it fails the binding criteria

Criteria Details
****************

.. list-table::

   * - Binding Criteria
     - Pass if Best Peptide is strong binder
     - ``IC50 MT < binding_threshold`` and ``%ile MT < percentile_threshold``
       (if ``--percentile-threshold`` parameter is set and 'conservative' ``--percentile-threshold-strategy`` is used) or
       ``IC50 MT < binding_threshold`` or ``%ile MT < percentile_threshold``
       (if 'exploratory' ``--percentile-threshold-strategy`` is used)
   * - Read Support Criteria
     - Pass if the variant has read support
     - ``Read Support < read_support``
   * - Expression Criteria
     - Pass if Best Transcript is expressed
     - ``Expr < expn_val``


.. _pvacfuse_reference_matches:

aggregated.tsv.reference_matches Report Columns
-----------------------------------------------

This file is only generated when the ``--run-reference-proteome-similarity``
option is chosen.

.. flat-table::
   :header-rows: 1

   * - Column Name
     - Description (BLAST)
     - Description (reference fasta)
   * - ``ID``
     - :cspan:`2` A unique identifier for the fusion
   * - ``Epitope Seq``
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

