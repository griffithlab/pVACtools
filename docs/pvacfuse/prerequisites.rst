.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Prerequisites
=============

Fusion detection and annotation
-------------------------------

pVACfuse expects an annotated bedpe file with fusion information. Fusion
detection will be preformed using `INTEGRATE <https://sourceforge.net/p/integrate-fusion/wiki/Home>`_ 
with annotations from `INTEGRATE-Neo <https://github.com/ChrisMaherLab/INTEGRATE-Neo>`_.

1. Align RNA with Tophat2 to obtain accepted_hits.bam and unmapped.bam
2. (OPTIONAL) Align WGS DNA with BWA aln/sampe (NOT MEM) to obtain tumor.dna.bam and normal.dna.bam
3. Produce gene annotations file with `gtfToGenePred <https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html>`_

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
    /gscmnt/gc2502/griffithlab/yafeng -e 8,9,10,11 fusions.bedpe.annot sample
    HLA-A*29:02,HLA-A*29:02,HLA-B*08:01,HLA-B*45:01,HLA-C*07:01,HLA-C*06:02
    NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC
    SMMalign hcc1395_fuse


.. Describe how to install and run INTEGRATE-Neo
.. Describe input file format
