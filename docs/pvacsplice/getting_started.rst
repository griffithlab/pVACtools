Getting Started
===============

pVACsplice provides a set of example data to show the expected format of input and output files. You can download the data set by running the ``pvacsplice download_example_data`` command.

The example data output can be reproduced by running the following command:

.. code-block:: none

   pvacsplice run \
   <example_data_dir>/HCC1395.splice_junctions.tsv \
   HCC1395_TUMOR_DNA \
   HLA-A*29:02,HLA-B*45:01,DRB1*04:05 \
   all \
   <output_dir> \
   <example_data_dir>/annotated.expression.vcf.gz \
   <example_data_dir>/ref_genome.fa \
   <example_data_dir>/Homo_sapiens.GRCh38.105.chr.gtf.gz \
   -e1 8,9,10 \
   -e2 15 \
   --normal-sample-name HCC1395_NORMAL_DNA

A detailed description of all command options can be found on the following page.
