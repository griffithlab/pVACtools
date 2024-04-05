.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Adding genotype sample information to your VCF
==============================================

pVACseq was primarily designed for clinical application. As such, it requires
that the input VCF contains sample genotype information (``GT`` field), which identifies
whether or not a variant was called in a specific sample of interest.

Some somatic variant callers (e.g., Strelka), however, do not include this field. In
other use cases you might want to run pVACseq on a list of variants of interest. If your
input VCF does not contain sample information (i.e. no ``FORMAT`` column
and/or sample column) or the ``FORMAT`` list does not contain a ``GT`` field,
you will need to preprocess your VCF to add this information.

This information can be added using the `VAtools <http://www.vatools.org>`_
``vcf-genotype-annotator``.

Using the vcf-genotype-annotator to add genotype information to your VCF
------------------------------------------------------------------------

Installing the vcf-genotype-annotator
*************************************

The ``vcf-genotype-annotator`` is part of the ``vatools`` package. 
Please visit `vatools.org <http://vatools.org>`_ for more details on this package.
You can install this package by running:

.. code-block:: none

   pip install vatools

Running the vcf-genotype-annotator
**********************************

**Example vcf-genotype-annotator commands**

.. code-block:: none

   vcf-genotype-annotator <input_vcf> <sample_name> 0/1 -o <gt_annotated_vcf>

The ``sample_name`` argument is used as the sample name in the ``#CHROM`` header
of your VCF when adding a new sample with this tool. If you want to add a ``GT``
field to an existing sample in your VCF, this argument will need to match the name of that sample.

Please see the `VAtools documentation
<https://vatools.readthedocs.io/en/latest/vcf_genotype_annotator.html>`_
for more information.
