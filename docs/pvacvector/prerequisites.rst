.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Prerequisites
=============

There are two options for the input file when running the pVACvector tool:

- A fasta file. This file contains protein sequences of candidate neoeptiopes
  to use for vector design.
- A :ref:`pvacseq` output TSV. This file has been filtered to include
  only the neoepitopes to use for vector design. If this file type is
  used, it is also necessary to provide the original VCF used in the
  pVACseq run via the ``--input_vcf`` option.
