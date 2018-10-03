.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Adding expression data to your VCF
==================================

pVACseq is able to parse coverage and expression information directly from the
VCF. The expected annotation format is outlined below.

===================== ==================================== =============================
Type                  VCF Sample                           Format Fields
===================== ==================================== =============================
Transcript Expression single-sample VCF or ``sample_name`` ``TX``
Gene Expression       single-sample VCF or ``sample_name`` ``GX``
===================== ==================================== =============================

**Transcript Expression**

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

**Gene Expression**

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

Using the vcf-expression-annotator to add expression information to your VCF
----------------------------------------------------------------------------

The ``vcf-expression-annotator`` will add expression information to your VCF.
It will accept expression data from various tools. Currently it supports
Cufflinks, Kallisto, StringTie, as well as a custom option for any
tab-delimited file.

Installing the vcf-expression-annotator
***************************************

The ``vcf-expression-annotator`` is part of the ``vcf-annotation-tools`` package.
You can install this package by running:

.. code-block:: none

   pip install vcf-annotation-tools

Running the vcf-expression-annotator
************************************

You can now use the output file from your expression caller to add expression information to
your VCF:

.. code-block:: none

    vcf-expression-annotator input_vcf expression_file kallisto|stringtie|cufflinks|custom gene|transcript

The data type ``gene`` or ``transcript`` identifies whether you are annotating
transcript or gene expression data. Transcript expression annotations will be
written to the ``TX`` format field while gene expression annotations will be
written to the ``GX`` format field. Please see the `documentation
<https://vcf-annotation-tools.readthedocs.io/en/latest/vcf_readcount_annotator.html>`_
for more information.
