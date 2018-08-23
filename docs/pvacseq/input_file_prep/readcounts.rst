.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Adding coverage data to your VCF
================================

pVACseq is able to parse coverage and information directly from the
VCF. The expected annotation format is outlined below.

===================== ==================================== =============================
Type                  VCF Sample                           Format Fields
===================== ==================================== =============================
Tumor DNA Coverage    single-sample VCF or ``sample_name`` ``AD``, ``DP``, and ``AF``
Tumor RNA Coverage    single-sample VCF or ``sample_name`` ``RAD``, ``RDP``, and ``RAF``
Normal DNA Coverage   ``--normal-sample-name``             ``AD``, ``DP``, and ``AF``
===================== ==================================== =============================

**Tumor DNA Coverage**

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the tumor DNA depth is determined from the ``DP`` format field.
The tumor DNA VAF is determined from the ``AF`` field. If the VCF does not contain a
``AF`` format field, the tumor DNA VAF is calculated from the ``AD`` and ``DP`` fields
by dividing the allele count by the total read depth.

**Tumor RNA Coverage**

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample the tumor RNA depth is determined from the ``RDP`` format field.
The tumor RNA VAF is determined from the ``RAF`` field. If the VCF does not contain a
``RAF`` format field, the tumor RNA VAF is calculated from the ``RAD`` and ``RDP`` fields
by dividing the allele count by the total read depth.

**Normal DNA Coverage**

To parse normal DNA coverage information, the input VCF to pVACseq will need to be a
multi-sample (tumor/normal) VCF, with one sample being the tumor sample, and the other
the matched normal sample. The tumor sample is identified by the
``sample_name`` parameter while the normal sample can be specified with
``--normal-sample-name`` option.

For this normal sample the normal DNA depth is determined from the ``DP`` format field.
The normal DNA VAF is determined from the ``AF`` field. If the VCF does not contain a
``AF`` format field, the normal DNA VAF is calculated from the ``AD`` and ``DP`` fields
by dividing the allele count by the total read depth.

Using the vcf-readcount-annotator to add coverage information to your VCF
-------------------------------------------------------------------------

Some variant callers will already have added coverage information to your VCF.
However, if your VCF doesn't contain coverage information or if you need to
add coverage information for additional samples or for RNA-seq data, you can
use the ``vcf-readcount-annotator`` to do so. The ``vcf-readcount-annotator``
will take the output from `bam-readcount
<https://github.com/genome/bam-readcount#build-instructions>`_ and use it to
add readcounts to your VCF.

Installing bam-readcount
************************

The ``vcf-readcount-annotator`` will add readcount information from bam-readcount
output files to your VCF. Therefore, you will first need to run bam-readcount
to obtain a file of readcounts for your variants.

Follow the installation instructions on the
`bam-readcount GitHub page <https://github.com/genome/bam-readcount#build-instructions>`_.

Running bam-readcount
*********************

bam-readcount uses a bam file and regions file as input, and the bam regions may
either contain snvs or indels. Indel regions must be run in a special insertion-centric
mode. Any mixed input regions must be split into snvs and indels, and bam-reacount must
then be run on each file individually using the same bam.

**Example bam-readcount command**

.. code-block:: none

   bam-readcount -f <reference fasta> -l <site list> <bam_file> [-i] [-b 20]

The ``-i`` option must be used when running indels bam in order to process indels in
insertion-centric mode.

A minimum base quality of 20 is recommended which can be enabled by ``-b 20``.

Installing the vcf-readcount-annotator
**************************************

The ``vcf-readcount-annotator`` is part of the ``vcf-annotation-tools`` package.
You can install this package by running:

.. code-block:: none

   pip install vcf-annotation-tools

Running the vcf-readcount-annotator
***********************************

You can now use the bam-readcount output file to add readcount information to
your VCF:

.. code-block:: none

   vcf-readcount-annotator input_vcf bam_readcount_file DNA|RNA -s sample_name

The data type ``DNA`` or ``RNA`` identifies whether you are annotating DNA or RNA
readcount. DNA readcount annotations will be written to the ``AD/DP/AF``
format fields while RNA readcount annotations will be written to the
``RAD/RDP/RAF`` format fields. Please see the `documentation
<https://vcf-annotation-tools.readthedocs.io/en/latest/vcf_readcount_annotator.html>`_
for more information.
