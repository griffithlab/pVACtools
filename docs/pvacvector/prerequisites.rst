.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Prerequisites
=============

There are two options for the input file when running the pVACvector tool:

- A FASTA file. This file contains protein sequences of candidate neoepitopes
  to use for vector design.
- A :ref:`pvacseq` output TSV. This file has been filtered to include
  only the neoepitopes to use for vector design. If this file type is
  used, it is also necessary to provide the original VCF used in the
  pVACseq run via the ``--input_vcf`` option. Output TSVs from MHC Class I and
  Class II pVACseq results can be combined into one by concatenating the two files and
  removing the duplicate header line.

Note that if you supply a FASTA file of peptides, these peptides will be used directly in the 
analysis and used in the final output. However, if you use a pvacseq TSV and variants VCF 
then the length of peptides extracted for junctional epitope testing and reporting in your output 
will be determined by the ``--input-n-mer`` option. 

