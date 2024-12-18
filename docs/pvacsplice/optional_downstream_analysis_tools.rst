.. image:: ../images/pVACsplice_logo_trans-bg_v4b.png
    :align: right
    :alt: pVACsplice logo
    :width: 175px

.. _pvacsplice_optional_downstream_analysis_tools_label:

Optional Downstream Analysis Tools
==================================

Generate Protein Fasta
----------------------

.. program-output:: pvacsplice generate_protein_fasta -h

This tool will extract protein sequences surrounding splice sites predicted by RegTools.
One use case for this tool is to help select long peptides that contain short neoepitope 
candidates. For example, if pVACsplice was run to predict nonamers (9-mers) that are good binders and
the user wishes to select long peptide (e.g. 24-mer) sequences that contain the nonamer for synthesis
or encoding in a DNA vector. The splice site junction will be centered in the
protein sequence returned (if possible).

The output may be limited to PASS variants only by setting the ``--pass`` only
flag. Additionally, variants can be limited to specific transcript biotypes
using the ``--biotypes`` parameters, which is set to only include ``protein_coding``
transcripts by default.

The output can be further limited to only certain variants by providing
a pVACsplice report file to the ``--input-tsv`` argument. Only the peptide sequences for the epitopes in the TSV
will be used when creating the FASTA. If this argument is an aggregated TSV
file, use the ``--aggregate-report-evaluation`` parameter to only include
peptide sequences for epitopes matching the chosen Evaluation(s).

Generate Aggregated Report
--------------------------

.. program-output:: pvacsplice generate_aggregated_report -h

This tool produces an aggregated version of the all_epitopes TSV. It finds the best-scoring
epitope for each splice site variant, and outputs additional binding affinity, expression, and
coverage information for that epitope. It also gives information about the
total number of well-scoring epitopes for each variant, the number of
transcripts covered by those epitopes, as well as the HLA alleles that those
epitopes are well-binding to. Lastly, the report will bin variants into tiers
that offer suggestions as to the suitability of variants for use in vaccines.
For a full definition of these tiers, see the pVACsplice :ref:`output file documentation <pvacsplice_aggregated>`.

Calculate Reference Proteome Similarity
---------------------------------------

.. program-output:: pvacsplice calculate_reference_proteome_similarity -h

This tool will find matches of the epitope candidates in the reference proteome and return the results in an output
TSV & reference_match file pair. It requires the input of a pVACplice run's fasta file in order to look up the larger
peptide sequence the epitope was derived from. Any substring of that peptide
sequence that matches against the reference proteome and is at least as long as the specified match length, will be
considered a hit. This tool also requires the user to provide a filtered.tsv,
all_epitopes.tsv or aggregated.tsv pVACsplice report file as an input and any
candidates in this input file will be searched for.

This tool may be either run with BLASTp using either the ``refseq_select_prot`` or ``refseq_protein`` database.
By default this option uses the BLAST API but users may :ref:`independently install BLASTp <blast>`. Alternatively, users
may provide a reference proteome fasta file and this tool will string match on
the entries of this fasta file directly. This approach is recommended, because
it is significantly faster than BLASTp. Reference proteome fasta files may be
downloaded from Ensembl. For example, the latest reference proteome fasta for human
can be downloaded from `this
link <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz>`_.

For more details on the generated reference_match file,
see the pVACsplice :ref:`output file documentation <pvacsplice_reference_matches>`.

NetChop Predict Cleavage Sites
------------------------------

.. program-output:: pvacsplice net_chop -h

This tool uses NetChop to predict cleavage sites for neoepitopes from a pVACsplice run's filtered/all_epitopes
TSV.  In its output, it adds to the TSV 3 columns: Best Cleavage Position, Best Cleavage Score, and a
Cleavage Sites list.  Typically this step is done in the pVACsplice run pipeline for the filtered output TSV
when specified. This tool provides a way to manually run this on pVACseq's generated filtered/all_epitopes
TSV files so that you can add this information when not present, if desired.

You can view more information about these columns for pVACsplice in the :ref:`output file documentation <all_ep_and_filtered>`.

NetMHCStab Predict Stability
----------------------------

.. program-output:: pvacsplice netmhc_stab -h

This tool uses NetMHCstabpan to add stability predictions for neoepitopes from a pVACsplice run's
filtered/all_epitopes TSV.  In its output, it adds to the TSV 4 columns: Predicted Stability, Half Life,
Stability Rank, and NetMHCStab Allele.  Typically this step is done in the pVACsplice run pipeline for the
filtered output TSV when specified.  This tool provides a way to manually run this on pVACseq's generated
filtered/all_epitopes TSV files so that you can add this information when not present if desired.

You can view more informatnion about these columns for pVACsplice in
the :ref:`output file documentation <pvacsplice_all_ep_and_filtered>`.

Identify Problematic Amino Acids
--------------------------------

.. program-output:: pvacsplice identify_problematic_amino_acids -h

This tool is used to identify positions in an epitope with an amino acid that
is problematic for downstream processing, e.g. vaccine manufacturing. Since
this can differ from case to case, this tool requires the user to specify which
amino acid(s) to consider problematic. This can be specified in one of three
formats:

.. list-table::

 * - ``amino_acid(s)``
   - One or more one-letter amino acid codes. Any occurrence of this amino acid string,
     regardless of the position in the epitope, is problematic. When specifying more than
     one amino acid, they will need to occur together in the specified order.
 * - ``amino_acid:position``
   - A one letter amino acid code, followed by a colon separator, followed by a positive
     integer position (one-based). The occurrence of this amino acid at the position
     specified is problematic., E.g. G:2 would check for a Glycine at the second position
     of the epitope. The N-terminus is defined as position 1.
 * - ``amino_acid:-position``
   - A one letter amino acid code, followed by a colon separator, followed by a negative
     integer position. The occurrence of this amino acid at the specified position from
     the end of the epitope is problematic. E.g., G:-3 would check for a Glycine at the
     third position from the end of the epitope. The C-terminus is defined as position -1.

You may specify any number of these problematic amino acid(s), in any
combination, by providing them as a comma-separated list.

This tool may be used with any filtered.tsv or all_epitopes.tsv pVACsplice report
file.
