.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Adding coverage data to your VCF
================================

pVACseq is able to parse coverage information directly from the
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

For this tumor sample, the tumor DNA depth is determined from the ``DP`` format field.
The tumor DNA VAF is determined from the ``AF`` field. If the VCF does not contain a
``AF`` format field, the tumor DNA VAF is calculated from the ``AD`` and ``DP`` fields
by dividing the allele count by the total read depth.

**Tumor RNA Coverage**

If the VCF is a single-sample VCF, pVACseq assumes that this sample is the
tumor sample. If the VCF is a multi-sample VCF, pVACseq will look for the
sample using the ``sample_name`` parameter and treat that sample as the tumor
sample.

For this tumor sample, the tumor RNA depth is determined from the ``RDP`` format field.
The tumor RNA VAF is determined from the ``RAF`` field. If the VCF does not contain a
``RAF`` format field, the tumor RNA VAF is calculated from the ``RAD`` and ``RDP`` fields
by dividing the allele count by the total read depth.

**Normal DNA Coverage**

To parse normal DNA coverage information, the input VCF to pVACseq will need to be a
multi-sample (tumor/normal) VCF, with one sample being the tumor sample, and the other
the matched normal sample. The tumor sample is identified by the
``sample_name`` parameter while the normal sample can be specified with
``--normal-sample-name`` option.

For this normal sample, the normal DNA depth is determined from the ``DP`` format field.
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

bam-readcount needs to be run separately for snvs and indels so it is
recommended to first split multi-allelic sites by using a tool such as ``vt
decompose``.

Installing vt
*************

The ``vt`` tool suite can be installed by following the instructions `on their
page <https://genome.sph.umich.edu/wiki/Vt#Installation>`_.

Installing bam-readcount
************************

The ``vcf-readcount-annotator`` will add readcount information from bam-readcount
output files to your VCF. Therefore, you will first need to run bam-readcount
to obtain a file of readcounts for your variants.

Follow the installation instructions on the
`bam-readcount GitHub page <https://github.com/genome/bam-readcount#build-instructions>`_.

Installing the vcf-readcount-annotator
**************************************

The ``vcf-readcount-annotator`` is part of the ``vatools`` package. 
Please visit `vatools.org <http://vatools.org>`_ for more details on this package.
You can install this package by running:

.. code-block:: none

   pip install vatools

Running vt decompose
********************

**Example vt decompose command**

.. code-block:: none

   vt decompose -s <input_vcf> -o <decomposed_vcf>

Running bam-readcount
*********************

bam-readcount uses a bam file and site list regions file as input. The site lists are
created from your decomposed VCF, one for snvs and one for indels. Snvs and
indels are then run separately through bam-readcount using the same bam. Indel regions
must be run in a special insertion-centric mode.

**Example bam-readcount command**

.. code-block:: none

   bam-readcount -f <reference_fasta> -l <site_list> <bam_file> [-i] [-b 20]

The ``-i`` option must be used when running the indels site list in order to process indels in
insertion-centric mode.

A minimum base quality of 20 is recommended which can be enabled using the ``-b 20``
option.

The ``mgibio/bam_readcount_helper-cwl`` Docker container contains a
``bam_readcount_helper.py`` script that will create the snv and indel site list files
from a VCF and run bam-readcount. Information on that Docker container can be found here:
`dockerhub mgibio/bam_readcount_helper-cwl <https://hub.docker.com/r/mgibio/bam_readcount_helper-cwl>`_.

**Example bam_readcount_helper.py command**

.. code-block:: none

   /usr/bin/python /usr/bin/bam_readcount_helper.py \
   <decomposed_vcf> <sample_name> <reference_fasta> <bam_file> NOPREFIX <output_dir>

This will write two bam-readcount files to the ``<output_dir>``:
``<sample_name>_bam_readcount_snv.tsv`` and
``<sample_name>_bam_readcount_indel.tsv``, containing readcounts for the snv
and indel positions, respectively.

Running the vcf-readcount-annotator
***********************************

The readcounts for snvs and indels are then added to your VCF separately, by
running the ``vcf-readcount-annotator`` twice.

**Example vcf-readcount-annotator commands**

.. code-block:: none

   vcf-readcount-annotator <decomposed_vcf> <snv_bam_readcount_file> <DNA|RNA> \
   -s <sample_name> -t snv -o <snv_annotated_vcf>

   vcf-readcount-annotator <snv_annotated_vcf> <indel_bam_readcount_file> <DNA|RNA> \
   -s <sample_name> -t indel -o <annotated_vcf>

The data type ``DNA`` or ``RNA`` identifies whether you are annotating DNA or RNA
readcount. DNA readcount annotations will be written to the ``AD/DP/AF``
format fields while RNA readcount annotations will be written to the
``RAD/RDP/RAF`` format fields. Please see the `VAtools documentation
<https://vatools.readthedocs.io/en/latest/vcf_readcount_annotator.html>`_
for more information.
