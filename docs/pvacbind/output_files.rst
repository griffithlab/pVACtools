.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

Output Files
============

The pVACbind pipeline will write its results in separate folders depending on
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
   * - ``<sample_name>.<MHC_I|MHC_II|Combined>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``.
   * - ``<sample_name>.<MHC_I|MHC_II|Combined>.filtered.tsv``
     - The above file after applying all filters, with cleavage site and stability
       predictions added.
   * - ``<sample_name>.<MHC_I|MHC_II|Combined>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format. Not generated when only running with presentation or immunogenicity algorithms.
   * - ``<sample_name>.<MHC_I|MHC_II|Combined>.all_epitopes.aggregated.tsv.reference_matches`` (optional)
     - A file outlining details of reference proteome matches

Filters applied to the filtered.tsv file
----------------------------------------

The filtered.tsv file is the all_epitopes file with the following filters
applied (in order):

- Binding Filter
- Top Score Filter

Please see the :ref:`Standalone Filter Commands<pvacbind_filter_commands>`
documentation for more information on each individual filter. The standalone
filter commands may be useful to reproduce the filtering or to chose different
filtering thresholds.

.. _pvacbind_all_ep_and_filtered:

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Index``
     - The FASTA ID of the peptide sequence the epitope belongs to
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Sub-peptide Position``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Epitope Seq``
     - The epitope sequence
   * - ``Median IC50 Score``
     - Median IC50 binding affinity of the epitope of all prediction algorithms used
   * - ``Best IC50 Score``
     - Lowest IC50 binding affinity of all prediction algorithms used
   * - ``Best IC50 Score Method``
     - Prediction algorithm with the lowest IC50 binding affinity for this epitope
   * - ``Median Percentile``
     - Median percentile rank of the epitope of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile``
     - Lowest percentile rank of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile Method``
     - Prediction algorithm with the lowest percentile rank for this epitope
   * - ``Median IC50 Percentile``
     - Median binding percentile rank of the epitope of all binding prediction algorithms used (those that provide percentile output)
   * - ``Best IC50 Percentile``
     - Lowest binding percentile rank of all binding prediction algorithms used (those that provide percentile output)
   * - ``Best IC50 Percentile Method``
     - Binding prediction algorithm with the lowest binding percentile rank for this epitope
   * - ``Median Immunogenicity Percentile``
     - Median immunogenicity percentile rank of the epitope of all
       immunogenicity prediction algorithms used (those that provide percentile output)
   * - ``Best Immunogenicity Percentile``
     - Lowest immunogenicity percentile rank of all immunogenicity prediction algorithms used (those that provide percentile output)
   * - ``Best Immunogenicity Percentile Method``
     - Immunogenicity prediction algorithm with the lowest immunogenicity percentile rank for this epitope
   * - ``Median Presentation Percentile``
     - Median presentation percentile rank of the epitope of all presentatio prediction algorithms used (those that provide percentile output)
   * - ``Best Presentation Percentile``
     - Lowest presentation percentile rank of all presentatio prediction algorithms used (those that provide percentile output)
   * - ``Best Presentation Percentile Method``
     - Presentation prediction algorithm with the lowest presentation percentile rank for this epitope
   * - ``Individual Prediction Algorithm Scores and Percentiles`` (multiple)
     - IC50 binding affinity scores, binding scores, presentation scores, processing scores, or immunogenicity scores as well as percentile ranks
       for the ``Epitope Seq`` for the individual prediction algorithms used. Percentile scores may be ``NA`` if not
       provided by the prediction algorithm.
   * - ``Problematic Positions`` (optional)
     - A list of positions in the ``Epitope Seq`` that match the problematic amino
       acids defined by the ``--problematic-amino-acids`` parameter
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

.. _pvacbind_aggregated:

all_epitopes.aggregated.tsv Report Columns
--------------------------------------------

The ``all_epitopes.aggregated.tsv`` file is an aggregated version of the all_epitopes TSV.
It shows the best-scoring epitope
for each variant, and outputs binding affinity and other information for that epitope. It gives information about the
total number of well-scoring epitopes for each variant as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.

Only epitopes meeting the ``--aggregate-inclusion-binding-threshold`` are included in this report (default: 5000).
If the number of unique epitopes for a mutation meeting this threshold exceeds the
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
     - A unique identifier for the variant
   * - HLA Alleles (multiple)
     - For each HLA allele in the run, the number of this variant's epitopes that bound well
       to the HLA allele (with median binding affinity < 1000)
   * - ``Best Peptide``
     - The best epitope sequence (see Best Peptide Criteria
       below for more details on how this is determined)
   * - ``Allele``
     - The Allele that the Best Peptide is binding to
   * - ``Pos``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Prob Pos``
     - A list of positions in the Best Peptide that are problematic. ``None`` if the ``-–problematic-pos`` parameter was not set during the pVACfuse run
   * - ``Num Included Peptides``
     - The number of included peptides according to the
       ``--aggregate-inclusion-binding-threshold`` and
       ``--aggregate-inclusion-count-limit``
   * - ``Num Passing Peptides``
     - The number of included peptides for this mutation that are well-binding.
   * - ``IC50 MT``
     - Median or lowest IC50 binding affinity of the Best Peptide across all prediction algorithms used
   * - ``%ile MT``
     - Median or lowest percentile rank of the Best Peptide across all prediction algorithms used
   * - ``IC50 %ile MT``
     - Median or lowest binding percentile rank of the Best Peptide across all binding prediction algorithms used
   * - ``IM %ile MT``
     - Median or lowest immunogenicity percentile rank of the Best Peptide across all immunogenicity prediction algorithms used
   * - ``Pres %ile MT``
     - Median or lowest presentation percentile rank of the Best Peptide across all presentation prediction algorithms used
   * - ``Ref Match`` (T/F) (optional)
     - Was there a match of the peptide sequence to the reference proteome?
   * - ``Tier``
     - A tier suggesting the suitability of variants for use in vaccines.
   * - ``Evaluation``
     - Column to store the evaluation of each variant. Either ``Accept``, ``Reject``, or ``Review``.

Best Peptide Criteria
_____________________

To determine the Best Peptide, all peptides meeting the
``--aggregate-inclusion-threshold`` and ``--aggregate-inclusion-count-limit``
(see above) for a variant are evaluated as follows:

- Pick the entries with no ``Problematic Positions``.
- For the remaining entries, calculate a rank for all the metrics specified
  via the ``--top-score-metric2`` parameter and sum them. Whether the lowest or median value
  is considered for each metric is controlled by the ``--top-score-metric`` parameter.
  Sort the remaining entries on this sum rank followed by the rank of the first
  ``--top-score-metric2`` specified (to break
  any ties in the sum rank). Select the highest sorted entry.

The pVACbind Aggregate Report Tiers
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
   * - ``--binding-percentile-threshold``
     - Use this threshold to filter epitopes on the IC50 %ile MT score.
     - 2.0
   * - ``--presentation-percentile-threshold``
     - Use this threshold to filter epitopes on the Pres %ile MT score.
     - 2.0
   * - ``--immunogenicity-percentile-threshold``
     - Use this threshold to filter epitopes on the IM %ile MT score.
     - 2.0
   * - ``--percentile-threshold-strategy``
     - Specify the candidate inclusion strategy. The ``conservative`` option requires a candidate to pass the
       binding threshold, the binding percentile threshold, the presentation percentile threshold, AND the
       immunogenicity percentile threshold. The ``exploratory`` option requires a candidate to pass EITHER the
       binding threshold, the binding percentile threshold, the presentation percentile threshold, OR the
       immunogenicity percentile threshold.
     - conservative
   * - ``--run-reference-proteome-similarity``
     - Set this flag in order to run reference proteome similarity analysis
       and enable ``RefMatch`` tiering. Use ``--blastp-path``, ``--blastp-db``,
       and ``--peptide-fasta`` parameters to configure your run.
     - False
   * - ``--problematic-amino-acids``
     - Configure this parameter in order to define amino acids problematic for
       the desired therapy delivery platform and enable ``ProbPos`` tiering.
     - None

Tiers
*****

Given the thresholds provided above, the Best Peptide is evaluated and binned
into tiers as follows:

.. list-table::
   :header-rows: 1

   * - Tier
     - Criteria
   * - ``Pass``
     - Best Peptide passes the scores, reference match, and problematic
       position criteria
   * - ``PoorBinder``
     - Best Peptide fails the binding criteria but passes the presentation, immunogenicity, reference match, and problematic position criteria
   * - ``PoorImmunogenicity``
     - Best Peptide fails the immunogenicity criteria but passes the binding, presentation, reference match, and problematic position criteria
   * - ``PoorPresentation``
     - Best Peptide fails the presentation criteria but passes the binding, immunogenicity, reference match, and problematic position criteria
   * - ``RefMatch``
     - Best Peptide fails the reference match criteria but passes the scores and problematic
       position criteria
   * - ``ProbPros``
     - Best Peptide fails the problematic position criteria but passes the scores and reference match criteria
   * - ``Poor``
     - Best Peptide doesn't fit in any of the above tiers, usually if it fails
       two or more criteria

Criteria Details
****************

.. list-table::
   :header-rows: 1

   * - Criteria
     - Description
     - Evaluation Logic
   * - Binding Criteria
     - Pass if Best Peptide is strong binder
     - binding score criteria: ``IC50 MT < binding_threshold``

       binding percentile score criteria: ``IC50 %ile MT < binding_percentile_threshold``

       ``conservative`` ``--percentile-threshold-strategy``: needs to pass
       BOTH the binding score criteria AND the binding percentile score criteria

       ``exploratory`` ``--percentile-threshold-strategy``: needs to pass
       EITHER the binding score criteria OR the binding percentile score criteria
   * - Presentation Criteria
     - Pass if the Best Peptide is presented by the MHC
     - ``Pres %ile MT < presentation_percentile_threshold``
   * - Immunogenicity Criteria
     - Pass if the Best Peptide is immunogenic
     - ``IM %ile MT < immunogenicity_percentile_threshold``
   * - Scores Criteria
     - Pass if the Best Peptide is a strong binder, presented by the MHC,
       and/or immunogenic
     - ``conservative`` ``--percentile-threshold-strategy``: needs to pass
       the binding criteria, the presentation criteria, AND the immunogenicity criteria

       ``exploratory`` ``--percentile-threshold-strategy``: needs to pass
       the binding criteria, the presentation criteria, OR the immunogenicity criteria
   * - Reference Match Criteria
     - Pass if there are no reference proteome matches
     - ``Ref Match == False``
   * - Problematic Position Criteria
     - Best Peptide does not contain a problematic amino acid as defined by the
       ``--problematic-amino-acids`` parameters
     - ``Prob Pos == None``

The pVACbind Aggregate Report Sorting
_____________________________________

The aggregate report is sorted as follows:

.. list-table::
   :header-rows: 1

   * - Sort Criteria
     - Sort Order
   * - ``Tier`` column
     - "Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation",
       "RefMatch", "ProbPos", "Poor"
   * - Sum of the ascending ranks of
       the metrics selected via the ``--top-score-metric2`` parameter (possible values:
       ``IC50 MT``, ``%ile MT``, ``IC50 %ile MT``, ``Pres %ile MT``; default: ``IC50 MT``,
       ``%ile MT``).
     - Ascending sum rank
   * - First metric specified in the ``--top-score-metric2`` as a tie breaker
       for identical sum ranks
     - Ascending rank
   * - ``ID`` column
     - Alphabetical

.. _pvacbind_reference_matches:

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
     - :cspan:`2` A unique identifier for the variant
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
