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
   * - ``Mutation``
     - The FASTA ID of the peptide sequence the epitope belongs to
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Sub-peptide Position``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Epitope Seq``
     - The epitope sequence
   * - ``Median Score``
     - Median ic50 binding affinity of the epitope of all prediction algorithms used
   * - ``Best Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Best Score Method``
     - Prediction algorithm with the lowest ic50 binding affinity for this epitope
   * - ``Median Percentile``
     - Median binding affinity percentile rank of the epitope of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile``
     - Lowest binding affinity percentile rank of all prediction algorithms used (those that provide percentile output)
   * - ``Best Percentile Method``
     - Prediction algorithm with the lowest binding affinity percentile rank for this epitope
   * - ``Individual Prediction Algorithm Scores and Percentiles`` (multiple)
     - ic50 binding affinity scores and percentiles for the ``Epitope Seq`` for the individual prediction algorithms used
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

.. _pvacbind_reference_matches:

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

.. _pvacbind_aggregated:

all_epitopes.aggregated.tsv Report Columns
--------------------------------------------

The ``all_epitopes.aggregated.tsv`` file is an aggregated version of the all_epitopes TSV.
Like the all_epitopes.tsv and filtered.tsv reports, in order to keep the outputs consistent,
pVACbind uses the same output columns as pVACseq for this file but some of the values will
be ``NA`` if a column doesn't apply to pVACbind.
This report presents the best-scoring (lowest binding affinity)
epitope for each variant and outputs additional binding affinity for that epitope.
It also gives information about the total number of well-scoring epitopes for each variant,
as well as the HLA alleles that those epitopes are well-binding to.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``ID``
     - A unique identifier for the variant
   * - ``HLA Alleles`` (multiple)
     - For each HLA allele in the run, the number of this variant's epitopes that bound well
       to the HLA allele (with median binding affinity < 1000)
   * - ``Gene``
     - ``NA``
   * - ``AA Change``
     - ``NA``
   * - ``Num Passing Transcripts``
     - ``NA``
   * - ``Best Peptide``
     - The best-binding epitope sequence (lowest median binding affinity)
   * - ``Pos``
     - ``NA``
   * - ``Num Passing Peptides``
     - The number of unique well-binding peptides for this mutation.
   * - ``IC50 MT``
     - Median IC50 binding affinity of the best-binding epitope across all prediction algorithms used
   * - ``IC50 WT``
     - ``NA``
   * - ``%ile MT``
     - Median binding affinity percentile rank of the best-binding epitope across all prediction algorithms used (those that provide percentile output)
   * - ``%ile WT``
     - ``NA``
   * - ``RNA Expr``
     - ``NA``
   * - ``RNA VAF``
     - ``NA``
   * - ``RNA Depth``
     - ``NA``
   * - ``DNA VAF``
     - ``NA``
   * - ``Tier``
     - ``NA``
   * - ``Evaluation``
     - Column to store the evaluation of each variant. Either ``Accept``, ``Reject``, or ``Review``.
