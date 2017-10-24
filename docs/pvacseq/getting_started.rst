Getting Started
---------------

pVACseq provides a set of example data to show the expected input and output files. You can download the data set by running the ``pvacseq download_example_data`` :ref:`command <example_data>`.

The example data output can be reproduced by running the following command:

.. code-block:: none

   pvacseq run \
   <example_data_dir>/input.vcf \
   Test \
   HLA-G*01:09,HLA-E*01:01,H2-IAb \
   NetMHC PickPocket NNalign <output_dir> \
   -e 9,10 \
   -i <example_data_dir>/additional_input_file_list.yaml --tdna-vaf 20 \
   --net-chop-method cterm --netmhc-stab \
   --top-score-metric=lowest -d full --keep-tmp-files

A detailed description of all command options can be found on the :ref:`Usage <pvacseq_run>` page.
