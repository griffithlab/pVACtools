.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _prerequisites_label:

Input File Preparation
======================

The main input file to the pVACseq pipeline is a VCF file. The VCF needs to
contain sample genotype information (``GT`` field). The VCF needs to be annotated 
with VEP to add transcript information.

If filtering on variant allele fractions (VAFs), depth, and expression values is 
desired, the VCF also needs to be annotated with this data.

Refer to the following sections for instructions on how to annotate your VCF with 
these data and how to produce a VCF for proximal variant analysis.

.. toctree::
   :maxdepth: 1
   :glob:

   input_file_prep/vep
   input_file_prep/readcounts
   input_file_prep/expression
   input_file_prep/proximal_vcf
