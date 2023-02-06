Input File Preparation
======================

The main input files to the pVACsplice pipeline are an annotated VCF file and a RegTools output file (tsv).

Step 1: Annotate VCF file
******************************

Please see the `pVACseq Input File Preparation <https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep.html>`_ for instructions on how to prepare an annotated VCF file.

Step 2: Run RegTools
********************

RegTools is a set of tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context. To run pVACsplice, you must first run RegTools ``cis-splice-effects identify`` to generate a list of non-canonical splicing junctions created from *cis*-acting regulatory variants. The output tsv file (-o option) is an input to pVACsplice. Here is an example command:

.. code-block:: none

   regtools cis-splice-effects identify \
   -o <output tsv file> \
   -s <rna_strand> \
   annotated_variants.vcf \ 
   rna_alignments.bam \ 
   ref.fa \
   annotations.gtf

Please see the `RegTools documentation
<https://regtools.readthedocs.io/en/latest/commands/cis-splice-effects-identify/>`_
for more information.



