Common Errors
-------------

Input VCF Sample Information
____________________________

**VCF contains more than one sample but sample_name is not set.**

pVACseq supports running with a multi-sample VCF as input. However, in this case it
requires the user to pick the sample to analyze, as only variants that are
called in the specified sample will be processed.

When running a multi-sample VCF the ``sample_name`` parameter is used to
identify which sample to analyze. Take, for example, the following ``#CHROM``
VCF header:

.. code-block:: none

   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR

This VCF contains two samples, ``NORMAL`` and ``TUMOR``. Use ``TUMOR`` as the
``sample_name`` parameter to process the tumor sample, and ``NORMAL`` to
process the normal sample.

If the input VCF only contains a single sample, the ``sample_name`` parameter
does not need to match the sample name in the VCF.

**sample_name not a sample ID in the #CHROM header of VCF**

This error occurs when running a multi-sample VCF and the ``sample_name``
parameter doesn't match any of the sample IDs in the VCF ``#CHROM`` header.
Take, for example, the following ``#CHROM`` header:

.. code-block:: none

   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR

All columns after ``FORMAT`` are sample identifiers that can be used as the
``sample_name`` parameter when running pVACseq, depending on which sample the
user wishes to process. Change the ``sample_name`` parameter of your ``pvacseq
run`` command to match one of them.

**normal_sample_name not a sample ID in the #CHROM header of VCF**

Your ``pvacseq run`` command included the ``--normal-sample-name`` parameter.
However, the argument chosen did not match any of the sample identifiers in
the ``#CHROM`` header of the input VCF.

Take, for example, the following ``#CHROM`` VCF header:

.. code-block:: none

   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR

All columns after ``FORMAT`` are sample identifiers that can be used as the
``--normal-sample-name`` parameter when running pVACseq, depending on which
sample is the normal sample in the VCF. Change the ``--normal-sample-name`` parameter of your ``pvacseq
run`` command to match the appropriate sample identifier.

**VCF doesn't contain any sample genotype information.**

pVACseq uses the sample genotype to identified which variants were called.
Therefore, while a VCF without a ``FORMAT`` and sample column(s) is valid, it cannot be used
in pVACseq. You will need to manually edit your VCF and add a ``FORMAT`` and
sample column with the ``GT`` genotype field. For more information on this
formatting please see the `VCF specification <https://github.com/samtools/hts-specs>`_ for your specific VCF version.

Input VCF Compression and Indexing
__________________________________

**Input VCF needs to be bgzipped when running with a proximal variants VCF.**

When running pVACseq with the ``--proximal-variants-vcf`` argument, the main
input VCF needs to be bgzipped and tabix indexed. See :ref:`the Input File
Preparation section <bgzip_input_vcf>` of the documentation for instructions on how to do so.

**Proximal variants VCF needs to be bgzipped.**

The VCF provided via the ``--proximal-variants-vcf`` argument needs to be
bgzipped and tabix indexed. See :ref:`the Input File
Preparation section <bgzip_phased_vcf>` of the documentation for instructions on how to do so.

**No .tbi file found for input VCF. Input VCF needs to be tabix indexed if processing with proximal variants.**

When running pVACseq with the ``--proximal-variants-vcf`` argument, the main
input VCF needs to be bgzipped and tabix indexed. See :ref:`the Input File
Preparation section <bgzip_input_vcf>` of the documentation for instructions on how to do so.

**No .tbi file found for proximal variants VCF. Proximal variants VCF needs to be tabix indexed.**

The VCF provided via the ``--proximal-variants-vcf`` argument needs to be
bgzipped and tabix indexed. See :ref:`the Input File
Preparation section <bgzip_phased_vcf>` of the documentation for instructions on how to do so.

Input VCF VEP Annotation
________________________

**Input VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.**

pVACseq requires the input VCF to be annotated by VEP. The provided input VCF
doesn't contain a ``CSQ`` ``INFO`` header. This indicates that it has not been
annotated. :ref:`The Input File Preparation section <vep>` of the
documentation provides instructions on how to annotate your VCF with VEP.

**VCF doesn't contain VEP DownstreamProtein annotations. Please re-annotate the VCF with VEP and the Wildtype and Downstream plugins.**

Although the input VCF was annotated with VEP, it is missing the required
annotations provided by the VEP Downstream plugin. The input VCF will need to
be reannotated using all of the required arguments as outlined in the :ref:`Input
File Preparation section <vep>` of the documentation.

**VCF doesn't contain VEP WildtypeProtein annotations. Please re-annotate the VCF with VEP and the Wildtype and Downstream plugins.**

Although the input VCF was annotated with VEP, it is missing the required
annotations provided by the VEP Wildtype plugin. The input VCF will need to
be reannotated using all of the required arguments as outlined in the :ref:`Input
File Preparation section <vep>` of the documentation.

**Proximal Variants VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.**

When running pVACseq with the ``--proximal-variants-vcf`` argument, that
proximal variants VCF needs to be annotated by VEP. The provided proximal
variants VCF
doesn't contain a ``CSQ`` ``INFO`` header. This indicates that it has not been
annotated. :ref:`The Input File Preparation section <vep>` of the
documentation provides instructions on how to annotate your VCF with VEP.

**There was a mismatch between the actual wildtype amino acid sequence and the expected amino acid sequence. Did you use the same reference build version for VEP that you used for creating the VCF?**

This error occurs when the reference nucleotide at a specific position is
different than the Ensembl transcript nucleotide at the same position. This results in
the mutant amino acid in the ``Amino_acids`` VEP annotation being different
from the amino acid of the transcript protein sequence as predicted by the
Wildtype plugin. The ``Amino_acids`` VEP annotation is based on the reference
and alternate nucleotides of the variant while the ``WildtypeProtein``
prediction is based on the Ensembl transcript nucleotide sequence.

This points to a fundamental disagreement between the reference that was
used during alignment and variant calling and the Ensembl reference. This
mismatch cannot be resolved by pVACseq, which is why this error is fatal.

Here are a few things that might resolve this error:

- If a VEP build 38 cache was accidentially used with a build 37 VCF (or vice
  versa), the correct cache needs to be downloaded and used during VEP annotation.
- Using the ``--assembly`` parameter during VEP annotation with the
  correct build version to match your VCF
- Using the ``fasta`` parameter during VEP annotation with the reference used
  to create the VCF
- Manually fixing the reference bases in your VCF to match the one expected by
  Ensembl
- Realigning and redoing variant calling on your sample with a reference that
  matches what is expected by VEP

If this mismatch cannot be resolved the VCF cannot be used by pVACseq.

Other
_____


**The TSV file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.**

None of the variants in the VCF file are supported by pVACseq.


.. A proximal variants TSV output path and peptide length need to be specified if a proximal variants input VCF is provided.")
   'Failed to extract format string from info description for tag (CSQ)')
   "Warning: TSV index already exists: {}".format(index))
   'Duplicate TSV indexes')
