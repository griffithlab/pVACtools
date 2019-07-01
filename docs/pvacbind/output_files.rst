.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

Output Files
============

The pVACbind pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms
- ``combined``: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoeptiope predictions from both

Each folder will contain the same list of output files (listed in the order
created):

=================================================== ===========
File Name                                           Description
=================================================== ===========
``<sample_name>.tsv``                               An intermediate file with variant and transcript information parsed from the input file.
``<sample_name>.tsv_<chunks>`` (multiple)           The above file but split into smaller chunks for easier processing with IEDB.
``<sample_name>.all_epitopes.tsv``                  A list of all predicted epitopes and their binding affinity scores, with additional variant information from the ``<sample_name>.tsv``.
``<sample_name>.filtered.tsv``                      The above file after applying all filters.
=================================================== ===========

Final Report Columns
--------------------

=============================================================== ===========
Column Name                                                     Description
=============================================================== ===========
``Mutation``                                                    The FASTA ID of the peptide sequence the epitope belongs to
``HLA Allele``                                                  The HLA allele for this prediction
``Sub-peptide Position``                                        The one-based position of the epitope in the protein sequence used to make the prediction
``Epitope Seq``                                                 The epitope sequence
``Median Score``                                                Median ic50 binding affinity of the epitope of all prediction algorithms used
``Best Score``                                                  Lowest ic50 binding affinity of all prediction algorithms used
``Best Score Method``                                           Prediction algorithm with the lowest ic50 binding affinity for this epitope
``Individual Prediction Algorithm  Scores`` (multiple)          ic50 scores for the ``Epitope Seq`` for the individual prediction algorithms used
``Best Cleavage Position`` (optional)                           Position of the highest predicted cleavage score
``Best Cleavage Score`` (optional)                              Highest predicted cleavage score
``Cleavage Sites`` (optional)                                   List of all cleavage positions and their cleavage score
``Predicted Stability Half Life`` (optional)                    The stability half life of the ``MT Epitope Seq``
``Stability Rank`` (optional)                                   The % rank stability of the ``MT Epitope Seq``
``NetMHCstab allele`` (optional)                                Nearest neighbor to the ``HLA Allele``. Used for NetMHCstab prediction
=============================================================== ===========
