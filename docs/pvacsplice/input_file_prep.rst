.. image:: ../images/pVACsplice_logo_trans-bg_v4b.png
    :align: right
    :alt: pVACsplice logo
    :width: 175px

.. _pvacsplice_prerequisites_label:

Input File Preparation
======================

The main input files to the pVACsplice pipeline are an annotated VCF file and a RegTools output file (tsv).

Step 1: Annotate VCF file
******************************

Please see the `pVACseq Input File Preparation <https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep.html>`_ for
instructions on how to prepare an annotated VCF file with the following differences:

- pVACsplice does not require annotations with the VEP Frameshift and Wildtype plugins so those parameters can be omitted when
  running VEP.
- The VEP parameter ``--transcript_version`` is required.
- pVACsplice does not use a phased VCF of proximal variants so that step can be skipped.

Step 2: Run RegTools
********************

`RegTools <https://regtools.readthedocs.io/en/latest/>`_ is a set of tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context. To run pVACsplice, you must first run RegTools ``cis-splice-effects identify`` to generate a list of non-canonical splicing junctions created from *cis*-acting regulatory variants. The output tsv file (-o option) is an input to pVACsplice. Here is an example command:

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



