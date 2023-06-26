.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _phased_vcf:

Creating a phased VCF of proximal variants
==========================================

By default, pVACseq will evaluate all somatic variants in the input VCF in
isolation. As a result, if a somatic variant of interest has other somatic
or germline variants in proximity, the calculated wildtype and mutant protein
sequences might be incorrect because the amino acid changes of those proximal
variants were not taken into account.

To solve this problem, we added a new option to pVACseq in the ``pvactools``
release 1.1. This option, ``--phased-proximal-variants-vcf``, can be
used to provide the path to a phased VCF of proximal variants in addition to
the normal input VCF. This VCF is then used to incorporate amino acid changes of nearby
variants that are in-phase to a somatic variant of interest. This results in
corrected mutant and wildtype protein sequences that account for proximal
variants.

At this time, this option only handles missense proximal variants but we are
working on a more comprehensive approach to this problem.

Note that if you do not perform the proximal variants step, you should manually 
review the sequence data for all candidates (e.g. in IGV) for proximal variants
and either account for these manually, or eliminate these candidates. Failure to 
do so may lead to inclusion of incorrect peptide sequences.

How to create the phased VCF of proximal variants
-------------------------------------------------

Input files
___________

- ``tumor.bam``: A BAM file of tumor reads
- ``somatic.vcf``: A VCF of somatic variants
- ``germline.vcf``: A VCF of germline variants
- ``reference.fa``: The reference FASTA file

Required tools
______________

- `Picard <https://broadinstitute.github.io/picard/>`_
- `GATK <https://software.broadinstitute.org/gatk/>`_
- `bgzip <http://www.htslib.org/doc/bgzip.html>`_
- `tabix <http://www.htslib.org/doc/tabix.html>`_

Create the reference dictionary
_______________________________

.. code-block:: none

   java -jar picard.jar CreateSequenceDictionary \
   R=reference.fa \
   O=reference.dict

Extract a tumor-only VCF from the somatic VCF
_____________________________________________

.. code-block:: none

    /usr/bin/java -Xmx4g -jar /opt/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R reference.fa
    --variant somatic.vcf
    --sample-name tumor_sample_name
    --output tumor_only.vcf

Index the tumor-only VCF
________________________

.. code-block:: none

   tabix -p vcf tumor_only.vcf


Update sample names
___________________

The sample names in the ``tumor.bam``, the newly-created ``tumor_only.vcf``, and the
``germline.vcf`` need to match. If they don't you need to edit the sample names
in the VCF files to match the tumor BAM file.

Combine somatic and germline variants using GATK’s CombineVariants
__________________________________________________________________

.. code-block:: none

    /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R reference.fa \
    --variant germline.vcf \
    --variant tumor_only.vcf \
    -o combined_somatic_plus_germline.vcf \
    --assumeIdenticalSamples

Sort combined VCF using Picard
______________________________

.. code-block:: none

   /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar SortVcf \
   I=combined_somatic_plus_germline.vcf \
   O=combined_somatic_plus_germline.sorted.vcf \
   SEQUENCE_DICTIONARY=reference.dict

Phase variants using GATK’s ReadBackedPhasing
_____________________________________________

Unfortunately, the tool used for this step is no longer available in current versions of GATK. We recommend using GATK 3.6.0 to run this step.

.. code-block:: none

   /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
   -T ReadBackedPhasing \
   -R reference.fa \
   -I tumor.bam \
   --variant combined_somatic_plus_germline.sorted.vcf \
   -L combined_somatic_plus_germline.sorted.vcf \
   -o phased.vcf

Annotate VCF with VEP
_____________________

.. code-block:: none

   ./vep \
   --input_file <phased.vcf> --output_file <phased.annotated.vcf> \
   --format vcf --vcf --symbol --terms SO --tsl\
   --hgvs --fasta <reference build FASTA file location> \
   --offline --cache [--dir_cache <VEP cache directory>] \
   --plugin Downstream --plugin Wildtype \
   [--dir_plugins <VEP_plugins directory>] [--pick] [--transcript_version]


.. _bgzip_phased_vcf:

bgzip and index the phased VCF
______________________________

.. code-block:: none

   bgzip -c phased.annotated.vcf > phased.annotated.vcf.gz
   tabix -p vcf phased.annotated.vcf.gz

The resulting ``phased.vcf.gz`` file can be used as the input to the
``--phased-proximal-variants-vcf`` option.

.. _bgzip_input_vcf:

bgzip and index the pVACseq main input VCF
__________________________________________

In order to use the ``--phased-proximal-variants-vcf`` option you will also
need to bgzip and index the VCF you plan on using as the main input VCF to
pVACseq. This is usually the same somatic.vcf used as input for creating the
phased proximal VCF after all the required and desired optional
preprocessing steps (e.g. VEP annotation, adding readcount and expression
data) have been done on that VCF. Bgzipping and indexing the fully
pre-processed somatic VCF is the last step.

.. code-block:: none

   bgzip -c input.vcf > input.vcf.gz
   tabix -p vcf input.vcf.gz
