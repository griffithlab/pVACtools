.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Output Files
============

The pVACseq pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms

Each folder will contain the same list of output files (listed in the order
created):

=================================================== ===========
File Name                                           Description
=================================================== ===========
``log``                                             A log file of the parameters used for this run
``<sample_name>.tsv``                               An intermediate file with variant, transcript, coverage, vaf, and expression information parsed from the input files.
``<sample_name>.tsv_<chunks>`` (multiple)           The above file but split into smaller chunks for easier processing with IEDB.
``<sample_name>.combined.parsed.tsv``               A modified version of the ``<sample_name>.tsv`` file with binding scores from IEDB added.
``<sample_name>.filtered.binding.tsv``              The above file after filtering by binding threshold.
``<sample_name>.filtered.coverage.tsv`` (optional)  The above file after filtering on coverage, VAF, and expression values (optional).
``<sample_name>.filtered.top.tsv`` (optional)       The above file after picking the top epitope for each variant (optional).
``<sample_name>.chop.tsv`` (optional)               The above file with cleavage site predictions added (optional).
``<sample_name>.stab.tsv`` (optional)               The above file with stability predictions added (optional).
``<sample_name>.final.tsv`` (optional)              The final output file after all filtering and optional steps.
=================================================== ===========

Final Report Columns
--------------------

=============================================================== ===========
Column Name                                                     Description
=============================================================== ===========
``Chromosome``                                                  The chromosome of this variant
``Start``                                                       The start position of this variant in the zero-based, half-open coordinate system
``Stop``                                                        The stop position of this variant in the zero-based, half-open coordinate system
``Reference``                                                   The reference allele
``Variant``                                                     The alt allele
``Transcript``                                                  The Ensembl ID of the affected transcript
``Ensembl Gene ID``                                             The Ensembl ID of the affected gene
``Variant Type``                                                The type of variant. ``missense`` for missense mutations, ``inframe_ins`` for inframe insertions, ``inframe_del`` for inframe deletions, and ``FS`` for frameshift variants
``Mutation``                                                    The amnio acid change of this mutation
``Protein Position``                                            The protein position of the mutation
``Gene Name``                                                   The Ensembl gene name of the affected gene
``HLA Allele``                                                  The HLA allele for this prediction
``Peptide Length``                                              The peptide length of the epitope
``Sub-peptide Position``                                        The one-based position of the epitope in the protein sequence used to make the prediction
``Mutation Position``                                           The one-based position of the start of the mutation in the epitope. ``0`` if the start of the mutation is before the epitope
``MT Epitope Seq``                                              Mutant epitope sequence
``WT Epitope Seq``                                              Wildtype (reference) epitope sequence at the same position in the full protein sequence. ``NA`` if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated
``Best MT Score Method``                                        Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
``Best MT Score``                                               Lowest ic50 binding affinity of all prediction algorithms used
``Corresponding WT Score``                                      ic50 binding affinity of the wildtype epitope. ``NA`` if there is no ``WT Epitope Seq``.
``Corresponding Fold Change``                                   ``Corresponding WT Score`` / ``Best MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
``Tumor DNA Depth``                                             Tumor DNA depth at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor DNA VAF``                                               Tumor DNA variant allele frequency at this position. ``NA`` if VCF entry does not contain tumor DNA readcount annotation.
``Tumor RNA Depth``                                             Tumor RNA depth at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Tumor RNA VAF``                                               Tumor RNA variant allele frequency at this position. ``NA`` if VCF entry does not contain tumor RNA readcount annotation.
``Normal DNA Depth``                                            Normal DNA depth at this position. ``NA`` if VCF entry does not contain normal DNA readcount annotation.
``Normal DNA VAF``                                              Normal DNA variant allele frequency at this position. ``NA`` if VCF entry does not contain normal DNA readcount annotation.
``Gene Expression``                                             Gene expression value at this position. ``NA`` if VCF entry does not contain gene expression annotation.
``Transcript Expression``                                       Transcript expression value at this position. ``NA`` if VCF entry does not contain transcript expression annotation.
``Median MT Score``                                             Median ic50 binding affinity of the mutant epitope of all prediction algorithms used
``Median WT Score``                                             Median ic50 binding affinity of the wildtype epitope of all prediction algorithms used. ``NA`` if there is no ``WT Epitope Seq``.
``Median Fold Change``                                          ``Median WT Score`` / ``Median MT Score``. ``NA`` if there is no ``WT Epitope Seq``.
``Individual Prediction Algorithm WT and MT Scores`` (multiple) ic50 scores for the ``MT Epitope Seq`` and ``WT Eptiope Seq`` for the individual prediction algorithms used
``Best Cleavage Position`` (optional)                           Position of the highest predicted cleavage score
``Best Cleavage Score`` (optional)                              Highest predicted cleavage score
``Cleavage Sites`` (optional)                                   List of all cleavage positions and their cleavage score
``Predicted Stability Half Life`` (optional)                    The stability half life of the ``MT Epitope Seq``
``Stability Rank`` (optional)                                   The % rank stability of the ``MT Epitope Seq``
``NetMHCstab allele`` (optional)                                Nearest neighbor to the ``HLA Allele``. Used for NetMHCstab prediction
=============================================================== ===========
