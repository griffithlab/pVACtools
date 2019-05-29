.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Output Files
============

The pVACseq pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms
- ``combined``: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoeptiope predictions from both

Each folder will contain the same list of output files (listed in the order
created):

=================================================== ===========
File Name                                           Description
=================================================== ===========
``<sample_name>.tsv``                               An intermediate file with variant, transcript, coverage, vaf, and expression information parsed from the input files.
``<sample_name>.tsv_<chunks>`` (multiple)           The above file but split into smaller chunks for easier processing with IEDB.
``<sample_name>.all_epitopes.tsv``                  A list of all predicted epitopes and their binding affinity scores, with additional variant information from the ``<sample_name>.tsv``.
``<sample_name>.filtered.tsv``                      The above file after applying all filters, with cleavage site and stability predictions added.
``<sample_name>.filtered.condensed.ranked.tsv``     A condensed version of the filtered TSV with only the most important columns remaining, with a priority score for each neoepitope candidate added.
=================================================== ===========

all_epitopes.tsv and filtered.tsv Report Columns
------------------------------------------------

=============================================================== ===========
Column Name                                                     Description
=============================================================== ===========
``Chromosome``                                                  The chromosome of this variant
``Start``                                                       The start position of this variant in the zero-based, half-open coordinate system
``Stop``                                                        The stop position of this variant in the zero-based, half-open coordinate system
``Reference``                                                   The reference allele
``Variant``                                                     The alt allele
``Transcript``                                                  The Ensembl ID of the affected transcript
``Transcript Support Level``                                    The `transcript support level (TSL) <https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl>`_ of the affected transcript. ``NA`` if the VCF entry doesn't contain TSL information.
``Ensembl Gene ID``                                             The Ensembl ID of the affected gene
``Variant Type``                                                The type of variant. ``missense`` for missense mutations, ``inframe_ins`` for inframe insertions, ``inframe_del`` for inframe deletions, and ``FS`` for frameshift variants
``Mutation``                                                    The amnio acid change of this mutation
``Protein Position``                                            The protein position of the mutation
``Gene Name``                                                   The Ensembl gene name of the affected gene
``HGVSc``                                                       The HGVS coding sequence variant name
``HGVSp``                                                       The HGVS protein sequence variant name
``HLA Allele``                                                  The HLA allele for this prediction
``Peptide Length``                                              The peptide length of the epitope
``Sub-peptide Position``                                        The one-based position of the epitope within the protein sequence used to make the prediction
``Mutation Position``                                           The one-based position of the start of the mutation within the epitope sequence. ``0`` if the start of the mutation is before the epitope
``MT Epitope Seq``                                              The mutant epitope sequence
``WT Epitope Seq``                                              The wildtype (reference) epitope sequence at the same position in the full protein sequence. ``NA`` if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated
``Best MT Score Method``                                        Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
``Best MT Score``                                               Lowest ic50 binding affinity of all prediction algorithms used
``Corresponding WT Score``                                      ic50 binding affinity of the wildtype epitope. ``NA`` if there is no ``WT Epitope Seq``.
``Corresponding Fold Change``                                   ``Corresponding WT Score`` / ``Best MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
``Tumor DNA Depth``                                             Tumor DNA depth at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor DNA VAF``                                               Tumor DNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor RNA Depth``                                             Tumor RNA depth at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Tumor RNA VAF``                                               Tumor RNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Normal DNA Depth``                                            Normal DNA depth at this position. ``NA`` if VCF entry does not contain normal DNA readcount annotation.
``Normal DNA VAF``                                              Normal DNA variant allele frequency (VAF) at this position. ``NA`` if VCF entry does not contain normal DNA readcount annotation.
``Gene Expression``                                             Gene expression value for the annotated gene containing the variant. ``NA`` if VCF entry does not contain gene expression annotation.
``Transcript Expression``                                       Transcript expression value for the annotated transcript containing the variant. ``NA`` if VCF entry does not contain transcript expression annotation.
``Median MT Score``                                             Median ic50 binding affinity of the mutant epitope across all prediction algorithms used
``Median WT Score``                                             Median ic50 binding affinity of the wildtype epitope across all prediction algorithms used. ``NA`` if there is no ``WT Epitope Seq``.
``Median Fold Change``                                          ``Median WT Score`` / ``Median MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
``Individual Prediction Algorithm WT and MT Scores`` (multiple) ic50 scores for the ``MT Epitope Seq`` and ``WT Eptiope Seq`` for the individual prediction algorithms used
``Best Cleavage Position`` (optional)                           Position of the highest predicted cleavage score
``Best Cleavage Score`` (optional)                              Highest predicted cleavage score
``Cleavage Sites`` (optional)                                   List of all cleavage positions and their cleavage score
``Predicted Stability Half Life`` (optional)                    The stability half life of the ``MT Epitope Seq``
``Stability Rank`` (optional)                                   The % rank stability of the ``MT Epitope Seq``
``NetMHCstab allele`` (optional)                                Nearest neighbor to the ``HLA Allele``. Used for NetMHCstab prediction
=============================================================== ===========

.. image:: ../images/output_file_columns.png
    :alt: pVACseq ouput file columns illustration

filtered.condensed.ranked.tsv Report Columns
--------------------------------------------

==================== ===========
Column Name          Description
==================== ===========
``Gene Name``        The Ensembl gene name of the affected gene.
``Mutation``         The amino acid change of this mutation.
``Protein Position`` The protein position of the mutation.
``HGVSc``            The HGVS coding sequence name.
``HGVSp``            The HGVS protein sequence name.
``HLA Allele``       The HLA allele for this prediction.
``MT Epitope Seq``   Mutant epitope sequence.
``MT IC50``          If ``--top-score-metric`` is set to ``lowest``, this corresponds to the ``Best MT Score`` in the full report. If ``--top-score-metric`` is set to ``median`` this corresponds to the ``Median MT Score`` in the full report.
``WT IC50``          If ``--top-score-metric`` is set to ``lowest``, this corresponds to the ``Corresponding WT Score`` in the full report. If ``--top-score-metric`` is set to ``median`` this corresponds to the ``Median WT Score`` in the full report.
``Fold Change``      If ``--top-score-metric`` is set to ``lowest``, this corresponds to the ``Corresponding Fold Change`` in the full report. If ``--top-score-metric`` is set to ``median`` this corresponds to the ``Median Fold Change`` in the full report.
``Tumor DNA Depth``  Tumor DNA depth at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor DNA VAF``    Tumor DNA variant allele frequency at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor RNA Depth``  Tumor RNA depth at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Tumor RNA VAF``    Tumor RNA variant allele frequency at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Gene Expression``  Gene expression value at this position. ``NA`` if VCF entry does not contain gene expression annotation.
``Rank``             A priority rank for the neoepitope (best = 1).
==================== ===========

.. _score:

The pVACseq Neoeptiope Priority Rank
____________________________________

Each of the following 4 criteria are assigned a rank-ordered value (worst = 1):

- B = ``MT IC50`` binding affinity, with the lowest being the best.
- F = ``Fold Change`` between MT and WT alleles, with the highest being the best.
- M = Mutant allele expression, calculated as (``Gene Expression`` * ``Tumor RNA VAF``), with the highest being the best.
- D = ``Tumor DNA VAF``, with the highest being the best.

A score is calculated from the above ranks with the following formula: ``B + F + (M * 2) + (D / 2)``. This score is then converted to a rank (best = 1).
