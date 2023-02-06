Getting Started
===============

pVACsplice provides a set of example data to show the expected format of input and output files. You can download the data set by running the ``pvacsplice download_example_data`` command.

The example data output can be reproduced by running the following command:

.. code-block:: none

   pvacsplice run \
   <example_data_dir>/input_junctions.tsv \
   <sample_name> \
   HLA-A*01:01,HLA-A*02:01,HLA-B*15:01,HLA-B*57:01,HLA-C*03:03,HLA-C*06:02 \
   all_class_i \
   <output_dir> \
   annotated_variants.vcf \
   ref.fa

A detailed description of all command options can be found on the following page.