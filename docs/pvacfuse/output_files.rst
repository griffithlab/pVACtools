.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Output Files
============

The pVACfuse pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms
- ``combined``: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoeptiope predictions from both

Each folder will contain the same list of output files (listed in the order
created):

.. list-table::
   :header-rows: 1

   * - File Name
     - Description
   * - ``<sample_name>.tsv``
     - An intermediate file with variant and transcript information parsed from the input file(s).
   * - ``<sample_name>.tsv_<chunks>`` (multiple)
     - The above file but split into smaller chunks for easier processing with IEDB.
   * - ``<sample_name>.all_epitopes.tsv``
     - A list of all predicted epitopes and their binding affinity scores, with
       additional variant information from the ``<sample_name>.tsv``.
   * - ``<sample_name>.filtered.tsv``
     - The above file after applying all filters, with cleavage site and stability
       predictions added.
   * - ``<sample_name>.filtered.condensed.ranked.tsv``
     - A condensed version of the filtered TSV with only the most important columns
       remaining, with a priority score for each neoepitope candidate added.

Filters applied to the filtered.tsv file
----------------------------------------

The filtered.tsv file is the all_epitopes file with the following filters
applied (in order):

- Binding Filter
- Top Score Filter

Please see the :ref:`Standalone Filter Commands<pvacfuse_filter_commands>`
documentation for more information on each individual filter. The standalone
filter commends may be useful to reproduce the filtering or to chose different
filtering thresholds.

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

In order to keep the outputs consistent, pVACfuse uses the same output columns
as pVACseq but some of the values will be ``NA`` if a column doesn't apply to
pVACfuse.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Chromosome``
     - The chromosome of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Start``
     - The start position of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Stop``
     - The stop position of the 5p and 3p portion of the fusion, separated by " / "
   * - ``Reference``
     - ``fusion``
   * - ``Variant``
     - ``fusion``
   * - ``Transcript``
     - The Ensembl IDs of the affected transcripts
   * - ``Transcript Support Level``
     - ``NA``
   * - ``Ensembl Gene ID``
     - ``NA``
   * - ``Variant Type``
     - The type of fusion. ``inframe_fusion`` for inframe fusions, ``frameshift_fusion`` for frameshift fusions
   * - ``Mutation``
     - ``NA``
   * - ``Protein Position``
     - The position of the fusion in the fusion protein sequence
   * - ``Gene Name``
     - The Ensembl gene names of the affected genes
   * - ``HGVSc``
     - ``NA``
   * - ``HGVSp``
     - ``NA``
   * - ``HLA Allele``
     - The HLA allele for this prediction
   * - ``Peptide Length``
     - The peptide length of the epitope
   * - ``Sub-peptide Position``
     - The one-based position of the epitope in the protein sequence used to make the prediction
   * - ``Mutation Position``
     - ``NA``
   * - ``MT Epitope Seq``
     - Mutant epitope sequence
   * - ``WT Epitope Seq``
     - ``NA``
   * - ``Best MT Score Method``
     - Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
   * - ``Best MT Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Corresponding WT Score``
     - ``NA``
   * - ``Corresponding Fold Change``
     - ``NA``
   * - ``Tumor DNA Depth``
     - ``NA``
   * - ``Tumor DNA VAF``
     - ``NA``
   * - ``Tumor RNA Depth``
     - ``NA``
   * - ``Tumor RNA VAF``
     - ``NA``
   * - ``Normal Depth``
     - ``NA``
   * - ``Normal VAF``
     - ``NA``
   * - ``Gene Expression``
     - ``NA``
   * - ``Transcript Expression``
     - ``NA``
   * - ``Median MT Score``
     - Median ic50 binding affinity of the mutant epitope of all prediction algorithms used
   * - ``Median WT Score``
     - ``NA``
   * - ``Median Fold Change``
     - ``NA``
   * - ``Individual Prediction Algorithm WT and MT Scores`` (multiple)
     - ic50 scores for the ``MT Epitope Seq`` and ``WT Epitope Seq`` for the individual prediction algorithms used
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

filtered.condensed.ranked.tsv Report Columns
--------------------------------------------

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - ``Gene Name``
     - The Ensembl gene names of the affected genes
   * - ``Mutation``
     - ``NA``
   * - ``Protein Position``
     - The position of the fusion in the fusion protein sequence
   * - ``HGVSc``
     - ``NA``
   * - ``HGVSp``
     - ``NA``
   * - ``HLA Allele``
     - The HLA allele for this prediction.
   * - ``Mutation Position``
     - ``NA``
   * - ``MT Epitope Seq``
     - Mutant epitope sequence.
   * - ``Median MT Score``
     - Median ic50 binding affinity of the mutant epitope across all prediction algorithms used
   * - ``Median WT Score``
     - ``NA``
   * - ``Median Fold Change``
     - ``NA``
   * - ``Best MT Score``
     - Lowest ic50 binding affinity of all prediction algorithms used
   * - ``Corresponding WT Score``
     - ``NA``
   * - ``Corresponding Fold Change``
     - ``NA``
   * - ``Tumor DNA Depth``
     - ``NA``
   * - ``Tumor DNA VAF``
     - ``NA``
   * - ``Tumor RNA Depth``
     - ``NA``
   * - ``Tumor RNA VAF``
     - ``NA``
   * - ``Gene Expression``
     - ``NA``
   * - ``Rank``
     - A priority rank for the neoepitope (best = 1).


The pVACfuse Neoeptiope Priority Rank
_____________________________________

The underlying formula for calculating the pVACfuse rank is the same as it is
for :ref:`rank`. However, since only the binding affinity is available for
fusion predictions, the pVACfuse simply ranks the neoeptiopes according to
their binding affinity, with the lowest being the best. If the ``--top-score-metric``
is set to ``median`` (default) the ``Median MT Score`` is used. If it
is set to ``lowest`` the ``Best MT Score`` is used.
