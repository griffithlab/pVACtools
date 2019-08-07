.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Optional Downstream Analysis Tools
==================================

Generate Protein Fasta
----------------------

.. program-output:: pvacfuse generate_protein_fasta -h

This tool will extract protein sequences surrounding fusion variant in an by parsing Integrate-Neo or AGFusion
output. One use case for this tool is to help select long peptides that contain short neoepitope 
candidates. For example, if pvacfuse was run to predict nonamers (9-mers) that are good binders and
the user wishes to select long peptide (e.g. 24-mer) sequences that contain the nonamer for synthesis
or encoding in a DNA vector. The fusion position will be centered in the protein sequence returned (if possible).
If the fusion causes a frameshift, the full downstream protein sequence will be returned unless the user specifies otherwise 
as described above.
