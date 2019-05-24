.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Prerequisites
=============

Fusion detection and annotation
-------------------------------

pVACfuse accepts two types of inputs, either an annotated bedpe file with
fusion information from `INTEGRATE-Neo <https://github.com/ChrisMaherLab/INTEGRATE-Neo>`_
or a output directory from `AGfusion <https://github.com/murphycj/AGFusion>`_ (recommended).

AGfusion
________

AGfusion allows a user to annotate output files from several fusion callers
using the ``agfusion batch`` command. The below example is for annotating the
output from the STAR-Fusion caller but many other fusion callers are supported.
For a full list see the `AGfusion documentation <https://github.com/murphycj/AGFusion#input-from-fusion-finding-algorithms>`_.

.. code-block:: none

   agfusion batch \
   -f <star_fusion_tsv> \
   -a starfusion \
   -db agfusion.homo_sapiens.87.db \
   - <output_directory> \
   --middlestart \
   --noncanonical

The ``--middlestar`` flag is required in order to use the ouput with pVACfuse.
This will indicate the fusion position in the fusion peptide sequence.

The ``--noncanonical`` flag is optional and can be used to annotate the fusion
with informations from all possible transcripts. By default only canonical
transcripts are used.

INTEGRATE-Neo
_____________

Fusion
detection will be preformed using `INTEGRATE <https://sourceforge.net/p/integrate-fusion/wiki/Home>`_ 
with annotations from `INTEGRATE-Neo <https://github.com/ChrisMaherLab/INTEGRATE-Neo>`_. It should be 
possible to start with fusions from another caller, convert the output to bedpe format, annotate the 
bedpe with INTEGRATE-Neo and then feed these candidates into pVACfuse.

1. Align RNA with Tophat2 (a requirement of INTEGRATE) to obtain accepted_hits.bam and unmapped.bam
2. (OPTIONAL) Align WGS DNA with BWA aln/sampe (NOT MEM, a requirement of INTEGRATE) to obtain tumor.dna.bam and normal.dna.bam
3. Produce a gene annotations file with `gtfToGenePred <https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html>`_

.. code-block:: none

    gtfToGenePred -genePredExt -geneNameAsName2 ref.gtf ref.genePred
    cut -f 1-10,12 ref.genePred > tmp.txt
    echo -e "#GRCh37.ensGene.name\tGRCh37.ensGene.chrom\tGRCh37.ensGene.strand\tGRCh37.ensGene.txStart\tGRCh37.ensGene.txEnd\tGRCh37.ensGene.cdsStart\tGRCh37.ensGene.cdsEnd\tGRCh37.ensGene.exonCount\tGRCh37.ensGene.exonStarts\tGRCh37.ensGene.exonEnds\tGRCh37.ensemblToGeneName.value" > annot.txt
    cat tmp.txt >> annot.txt

4. Run `INTEGRATE <https://sourceforge.net/p/integrate-fusion/wiki/Home>`_ to obtain fusions.bedpe

.. code-block:: none

    Integrate fusion ref.fa annot.txt bwts accepted_hits.bam unmappeds.bam [tumor.dna.bam normal.dna.bam | tumor.dna.bam]

5. Run `INTEGRATE-Neo <https://github.com/ChrisMaherLab/INTEGRATE-Neo>`_ to obtain annotated fusions bedpe file

.. code-block:: none

    integrate-neo.py -t hla.optitype -f fusions.bedpe -r ref.fa -g ref.genePred -k

.. <===== pVACfuse =====>
    pvacfuse run --net-chop-method cterm --netmhc-stab --iedb-install-directory
    IEDB_INSTALL_DIRECTORY -e 8,9,10,11 fusions.bedpe.annot sample
    HLA-A*29:02,HLA-A*29:02,HLA-B*08:01,HLA-B*45:01,HLA-C*07:01,HLA-C*06:02
    NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC
    SMMalign output_dir

.. Describe how to install and run INTEGRATE-Neo
.. Describe input file format
