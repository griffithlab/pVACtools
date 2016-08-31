.. _prerequisites_label:

Prerequisites
=============

VEP
---

The input to the pVAC-Seq pipeline is a VEP annotated single-sample VCF. In addition to the standard VEP annotations, pVAC-Seq also requires the annotations provided by the Downstream and Wildtype VEP plugins.

To create a VCF for use with pVAC-Seq follow these steps:

1. Download and install the VEP command line tool following `these instructions <http://useast.ensembl.org/info/docs/tools/vep/script/index.html>`_.
2. Download the VEP_plugins from their `GitHub repository <https://github.com/Ensembl/VEP_plugins>`_.
3. :ref:`Copy the Wildtype plugin<install_vep_plugin_label>` provided with the pVAC-Seq package to the folder with the other VEP_plugins:

.. code-block:: none

   pvacseq install_vep_plugin

4. Run VEP on the input vcf with at least the following options:

.. code-block:: none

   --format vcf
   --vcf
   --symbol
   --plugin Downstream
   --plugin Wildtype
   --terms SO

The ``--dir_plugins <VEP_plugins directory>`` option may need to be set depending on where the VEP_plugins were installed to.

Additional VEP options that might be desired can be found `here <http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html>`_.

**Example VEP Command**

.. code-block:: none

   perl variant_effect_predictor.pl \
   --input_file <input VCF> --format vcf --output_file <output VCF> \
   --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype \
   [--dir_plugins <VEP_plugins directory>]

Optional Preprocessing
----------------------

Coverage and expression data can be added to the pVAC-Seq processing by providing bam-readcount and/or Cufflinks output files as additional input files. These additional input files must be provided as a yaml file in the following structure:

.. code-block:: none

    gene_expn_file: <genes.fpkm_tracking file from Cufflinks>
    transcript_expn_file: <isoforms.fpkm_tracking file from Cufflinks>
    normal_snvs_coverage_file: <bam-readcount output file for normal BAM and snvs>
    normal_indels_coverage_file: <bam-readcount output file for normal BAM and indels>
    tdna_snvs_coverage_file: <bam-readcount output file for tumor DNA BAM and snvs>
    tdna_indels_coverage_file: <bam-readcount output file for tumor DNA BAM and indels>
    trna_snvs_coverage_file: <bam-readcount output file for tumor RNA BAM and snvs>
    trna_indels_coverage_file: <bam-readcount output file for tumor RNA BAM and indels>

Any file in this list is optional and its entry can be omitted. If no additional files exist than this yaml file is optional and can be omitted from the list of pvacseq arguments.

bam-readcount
^^^^^^^^^^^^^

pVAC-Seq optionally accepts bam-readcount files as inputs which are used to add coverage information (depth and VAF) for downstream filtering. Depth and VAF are calculated from the read counts of the reference allele and alternate allele.

Follow the installation instructions on the `bam-readcount GitHub page <https://github.com/genome/bam-readcount#build-instructions>`_.

bam-readcount uses a bam file as input and the bam may either contain snvs or indels. Indel bams must be run in a special insertion-centric mode. Any mixed input bams must be split into snvs and indels and bam-reacount must then be run on each file individually.

**Example bam-readcount command**

.. code-block:: none

   bam-readcount -f <reference fasta> -l <site list> <bam_file>

The ``-i`` option must be used when running indels bam in order to process indels in insertion-centric mode.

A minimum base quality of 20 is recommended which can be enabled by ``-b 20``.

Cufflinks
^^^^^^^^^

Cufflinks files can be added as optional inputs to pVAC-Seq to extract gene and transcript expression data for downstream filtering.

Installation instructions for Cufflinks can be found on their `GitHub page <https://github.com/cole-trapnell-lab/cufflinks#install-quick-start>`_.

**Example Cufflinks command**

.. code-block:: none

   cufflinks <sam_file>
