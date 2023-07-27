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
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with cleavage site and stability
       predictions added.
   * - ``<sample_name>.all_epitopes.aggregated.tsv``
     - An aggregated version of the ``all_epitopes.tsv`` file that gives information about
       the best epitope for each mutation in an easy-to-read format.
   * - ``<sample_name>.all_epitopes.aggregated.tsv.reference_matches`` (optional)
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

.. _pvacbind_all_ep_and_filtered:

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Mutation``
     - The FASTA ID of the peptide sequence the epitope belongs to
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Sub-peptide Position``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Epitope Seq``
     - The epitope sequence
   * - ``Median IC50 Score``
     - Median ic50 binding affinity of the epitope of all prediction algorithms used
   * - ``Best IC50 Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Best IC50 Score Method``
     - Prediction algorithm with the lowest ic50 binding affinity for this epitope
   * - ``Median Percentile``
     - Median binding affinity percentile rank of the epitope of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile``
     - Lowest binding affinity percentile rank of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile Method``
     - Prediction algorithm with the lowest binding affinity percentile rank for this epitope
   * - ``Individual Prediction Algorithm Scores and Percentiles`` (multiple)
     - ic50 binding affinity scores and percentiles for the ``Epitope Seq`` for the individual prediction algorithms used
   * - ``MHCflurryEL WT and MT Processing Score and Presentation Score and Percentile`` (optional)
     - MHCflurry elution processing score and presentation score and percentiles
       for the ``MT Epitope Seq`` and ``WT Epitiope Seq`` if the run included
       MHCflurryEL as one of the prediction algorithms
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

Only epitopes meeting the ``--aggregate-inclusion-threshold`` are included in this report (default: 5000).
Whether the median or the lowest binding affinity metrics are output in the ``IC50 MT``,
``%ile MT``, and columns is controlled by the ``--top-score-metric`` parameter.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``ID``
     - A unique identifier for the variant
   * - ``HLA Alleles`` (multiple)
     - For each HLA allele in the run, the number of this variant's epitopes that bound well
       to the HLA allele (with median binding affinity < 1000)
   * - ``Best Peptide``
     - The best-binding epitope sequence (lowest median binding affinity)
   * - ``Prob Pos``
     - A list of positions in the Best Peptide that are problematic. ``None`` if the ``-–problematic-pos`` parameter was not set during the pVACfuse run
   * - ``Num Passing Peptides``
     - The number of unique well-binding peptides for this mutation.
   * - ``IC50 MT``
     - Median IC50 binding affinity of the best-binding epitope across all prediction algorithms used
   * - ``%ile MT``
     - Median binding affinity percentile rank of the best-binding epitope across all prediction algorithms used (those that provide percentile output)
   * - ``Ref Match`` (T/F) (optional)
     - Was there a match of the peptide sequence to the reference proteome?
   * - ``Evaluation``
     - Column to store the evaluation of each variant. Either ``Accept``, ``Reject``, or ``Review``.

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
   * - ``--percentile-threshold``
     - When set, use this threshold to filter epitopes on the %ile MT score in addition to having to meet the binding threshold.
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
     - Best Peptide passes the binding criteria
   * - ``Poor``
     - Best Peptide fails the binding criteria

Criteria Details
****************

.. list-table::

   * - Binding Criteria
     - Pass if Best Peptide is strong binder
     - ``IC50 MT < binding_threshold`` and ``%ile MT < percentile_threshold``
       (if ``--percentile-threshold`` parameter is set)


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
