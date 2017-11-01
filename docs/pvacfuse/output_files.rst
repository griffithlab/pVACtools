.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Output Files
============

The pVACfuse pipeline will write its results in separate folders depending on
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
``<sample_name>.filtered.top.tsv`` (optional)       The above file after picking the top epitope for each variant (optional).
``<sample_name>.chop.tsv`` (optional)               The above file with cleavage site predictions added (optional).
``<sample_name>.stab.tsv`` (optional)               The above file with stability predictions added (optional).
``<sample_name>.final.tsv`` (optional)              The final output file after all filtering and optional steps.
=================================================== ===========

Final Report Columns
--------------------

In order to keep the outputs consistent, pVACfuse uses the same output columns
as pVACseq but some of the values will be ``NA`` if a column doesn't apply to
pVACfuse.

=============================================================== ===========
Column Name                                                     Description
=============================================================== ===========
``Chromosome``                                                  The chromosome of the 5p and 3p portion of the fusion, separated by " / "
``Start``                                                       The start position of the 5p and 3p portion of the fusion, separated by " / "
``Stop``                                                        The stop position of the 5p and 3p portion of the fusion, separated by " / "
``Reference``                                                   ``fusion``
``Variant``                                                     ``fusion``
``Transcript``                                                  The Ensembl IDs of the affected transcripts
``Ensembl Gene ID``                                             ``NA``
``Variant Type``                                                The type fusion. ``inframe_fusion`` for inframe fusions, ``frameshift_fusion`` for frameshift fusions
``Mutation``                                                    ``NA``
``Protein Position``                                            The position of the fusion in the fusion protein sequence
``Gene Name``                                                   The Ensembl gene names of the affected genes
``HLA Allele``                                                  The HLA allele for this prediction
``Peptide Length``                                              The peptide length of the epitope
``Sub-peptide Position``                                        The one-based position of the epitope in the protein sequence used to make the prediction
``Mutation Position``                                           ``NA```
``MT Epitope Seq``                                              Mutant epitope sequence
``WT Epitope Seq``                                              ``NA``
``Best MT Score Method``                                        Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope
``Best MT Score``                                               Lowest ic50 binding affinity of all prediction algorithms used
``Corresponding WT Score``                                      ``NA``
``Corresponding Fold Change``                                   ``NA``
``Tumor DNA Depth``                                             ``NA``
``Tumor DNA VAF``                                               ``NA``
``Tumor DNA Depth``                                             ``NA``
``Tumor DNA VAF``                                               ``NA``
``Tumor DNA Depth``                                             ``NA``
``Tumor DNA VAF``                                               ``NA``
``Gene Expression``                                             ``NA``
``Transcript Expression``                                       ``NA``
``Median MT Score``                                             Median ic50 binding affinity of the mutant epitope of all prediction algorithms used
``Median WT Score``                                             ``NA``
``Median Fold Change``                                          ``NA``
``Individual Prediction Algorithm WT and MT Scores`` (multiple) ic50 scores for the ``MT Epitope Seq`` and ``WT Eptiope Seq`` for the individual prediction algorithms used
``Best Cleavage Position`` (optional)                           Position of the highest predicted cleavage score
``Best Cleavage Score`` (optional)                              Highest predicted cleavage score
``Cleavage Sites`` (optional)                                   List of all cleavage positions and their cleavage score
``Predicted Stability Half Life`` (optional)                    The stability half life of the ``MT Epitope Seq``
``Stability Rank`` (optional)                                   The % rank stability of the ``MT Epitope Seq``
``NetMHCstab allele`` (optional)                                Nearest neighbor to the ``HLA Allele``. Used for NetMHCstab prediction
=============================================================== ===========
