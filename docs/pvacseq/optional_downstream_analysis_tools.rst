.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _optional_downstream_analysis_tools_label:

Optional Downstream Analysis Tools
==================================

Generate Protein Fasta
----------------------

.. program-output:: pvacseq generate_protein_fasta -h

.. .. argparse::
    :module: lib.generate_protein_fasta
    :func: define_parser
    :prog: pvacseq generate_protein_fasta

This tool will extract protein sequences surrounding supported protein altering variants in an
input VCF file. One use case for this tool is to help select long peptides that contain short neoepitope 
candidates. For example, if pvacseq was run to predict nonamers (9-mers) that are good binders and
the user wishes to select long peptide (e.g. 24-mer) sequences that contain the nonamer for synthesis
or encoding in a DNA vector. The protein sequence extracted will correspond to the transcript sequence 
used in the annotated VCF. The alteration in the VCF (e.g. a somtic missense SNV) will be centered in the 
protein sequence returned (if possible). If the variant is near the beginning or end of the CDS, it will
be as close to center as possible while returning the desired protein sequence length. If the variant
causes a frameshift, the full downstream protein sequence will be returned unless the user specifies otherwise 
as described above. 

Generate Condensed, Ranked Report
---------------------------------

.. program-output:: pvacseq generate_condensed_ranked_report -h

This tool will produce a condensed version of the filtered TSV with only the most important columns remaining, 
with a score for each neoepitope candidate added. Refer to the Output Files section for more details on the 
format of this report.

