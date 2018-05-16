.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _prerequisites_label:

Prerequisites
=============

VEP
---

The input to the pVACseq pipeline is a VEP annotated VCF. In addition to the standard VEP annotations, pVACseq also requires the annotations provided by the Downstream and Wildtype VEP plugins.

To create a VCF for use with pVACseq follow these steps:

1. Download and install the VEP command line tool following `these instructions <http://useast.ensembl.org/info/docs/tools/vep/script/index.html>`_.
2. Download the VEP_plugins from their `GitHub repository <https://github.com/Ensembl/VEP_plugins>`_.
3. :ref:`Copy the Wildtype plugin<install_vep_plugin_label>` provided with the pVACseq package to the folder with the other VEP_plugins:

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

The ``--pick`` option might be useful to limit the annotation to the top
transcripts. Otherwise, VEP will annotate each variant with all possible
transcripts. pVACseq will provide predictions for all transcripts in the VEP
CSQ field. Running VEP without the ``--pick`` option can therefor drasticly
increase the runtime of pVACseq.

Additional VEP options that might be desired can be found `here <http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html>`_.

**Example VEP Command**

.. code-block:: none

   perl variant_effect_predictor.pl \
   --input_file <input VCF> --format vcf --output_file <output VCF> \
   --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype \
   [--dir_plugins <VEP_plugins directory>]

Optional Preprocessing to add coverage and expression data
----------------------------------------------------------

pVACseq is able to parse coverage and expression information directly from the
VCF. The expected annotation format is outlined below.

===================== ==================================== =============================
Type                  VCF Sample                           Format Fields
===================== ==================================== =============================
Tumor DNA Coverage    single-sample VCF or ``sample_name`` ``AD``, ``DP``, and ``AF``
Tumor RNA Coverage    single-sample VCF or ``sample_name`` ``RAD``, ``RDP``, and ``RAF``
Normal DNA Coverage   ``--normal-sample-name``             ``AD``, ``DP``, and ``AF``
Transcript Expression single-sample VCF or ``sample_name`` ``TX``
Gene Expression       single-sample VCF or ``sample_name`` ``GX``
===================== ==================================== =============================

Tumor DNA Coverage
^^^^^^^^^^^^^^^^^^

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the tumor DNA depth is determined from the ``DP`` format field.
The tumor DNA VAF is determined from the ``AF`` field. If the VCF does not contain a
``AF`` format field, the tumor DNA VAF is calculated from the ``AD`` and ``DP`` fields
by dividing the allele count by the total read depth.

Tumor RNA Coverage
^^^^^^^^^^^^^^^^^^

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the tumor RNA depth is determined from the ``RDP`` format field.
The tumor RNA VAF is determined from the ``RAF`` field. If the VCF does not contain a
``RAF`` format field, the tumor RNA VAF is calculated from the ``RAD`` and ``RDP`` fields
by dividing the allele count by the total read depth.

Normal DNA Coverage
^^^^^^^^^^^^^^^^^^^

To parse normal DNA coverage information, the input VCF to pVACseq will need to be a
multi-sample (tumor/normal) VCF, with one sample being the tumor sample, and the other
the matched normal sample. The tumor sample is identified by the
``sample_name`` parameter while the normal sample can be specified with
``--normal-sample-name`` option.

For this normal sample the normal DNA depth is determined from the ``DP`` format field.
The normal DNA VAF is determined from the ``AF`` field. If the VCF does not contain a
``AF`` format field, the normal DNA VAF is calculated from the ``AD`` and ``DP`` fields
by dividing the allele count by the total read depth.

Transcript Expression
^^^^^^^^^^^^^^^^^^^^^

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the transcript expression is determined from the ``TX``
format field. The ``TX`` format field is a comma-separated list of
per-transcript expression values, where each individual transcript expression
is listed as ``expression_id|expression_value``, e.g.
``ENST00000215794|2.35912,ENST00000215795|0.2``. The ``expression_id`` needs
to match the ``Feature`` field of the VEP ``CSQ`` annotation.

Gene Expression
^^^^^^^^^^^^^^^

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the gene expression is determined from the ``GX``
format field. The ``GX`` format field is a comma-separated list of
per-gene expression values, where each individual gene expression
is listed as ``gene_id|expression_value``, e.g.
``ENSG00000184979|2.35912``. The ``gene_id`` needs to match the ``Gene`` field
of the VEP ``CSQ`` annotation.

Adding coverage and expression annotations to your input VCF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We've developed the ``vcf-annotation-tools`` package that contains several
tools to add coverage and expression information to your VCFs. You can install
this package by running:

.. code-block:: none

   pip install vcf-annotation-tools

This command will install two tools, the ``vcf-readcount-annotator`` and the
``vcf-expression-annotator`` that can be used to add coverage annotations and
expression annotations to your VCF, respectively.


Running vcf-readcount-annotator to add coverage information to your input VCF
********************************************************************************

The ``vcf-readcount-annotator`` will add readcount information from
bam-readcount to your VCF.

Follow the installation instructions on the `bam-readcount GitHub page <https://github.com/genome/bam-readcount#build-instructions>`_.

bam-readcount uses a bam file and regions file as input, and the bam regions may either contain snvs or indels. Indel regions must be run in a special insertion-centric mode. Any mixed input regions must be split into snvs and indels, and bam-reacount must then be run on each file individually using the same bam.

**Example bam-readcount command**

.. code-block:: none

   bam-readcount -f <reference fasta> -l <site list> <bam_file>

The ``-i`` option must be used when running indels bam in order to process indels in insertion-centric mode.

A minimum base quality of 20 is recommended which can be enabled by ``-b 20``.

You can now use the bam-readcount output file to add readcount information to
your VCF:

.. code-block:: none

   vcf-readcount-annotator input_vcf bam_readcount_file DNA|RNA -s sample_name

The data type ``DNA`` or ``RNA`` identifies whether you are annotating DNA or RNA
readcount. DNA readcount annotations will be written to the ``AD/DP/AF``
format fields while RNA readcount annotations will be written to the
``RAD/RDP/RAF`` format fields.

Running vcf-expression-annotator to add expression information to your input VCF
********************************************************************************

The ``vcf-expression-annotator`` will add expression information to your VCF.
It will accept expression data from various tools. Currently it supports
Cufflinks, Kallisto, StringTie, as well as a custom option for any
tab-delimited file.

.. code-block:: none

    vcf-expression-annotator input_vcf expression_file kallisto|stringtie|cufflinks|custom gene|transcript

The data type ``gene`` or ``transcript`` identifies whether you are annotating
transcript or gene expression data. Transcript expression annotations will be
written to the ``TX`` format field while gene expression annotations will be
written to the ``GX`` format field.
