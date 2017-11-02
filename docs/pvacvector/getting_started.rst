.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Getting Started
---------------

pVACvector provides a set of example data to show the expected input and output files. You can download the data set by running the ``pvacfuse download_example_data`` :ref:`command <pvacvector_example_data>`.

There are two option as to how to run pVACvector depending on the input file
type used. You can either use a pVACseq output TSV of neoepitopes or a FASTA
file of peptide sequences.

Here is an example of how to run pVACvector with a pVACseq output TSV:

.. code-block:: none

   pvacvector run \
   <example_data_dir>/input.tsv \
   Test \
   H-2-Kb \
   NetMHC \
   <output_dir> \
   -e 8 \
   -v <example_data_dir>/input.vcf
   --keep-tmp-files

In this example pVACvector is run with an input FASTA file:

.. code-block:: none

   pvacvector run \
   <example_data_dir>/input.fa \
   Test \
   H-2-Kb \
   NetMHC \
   <output_dir> \
   -e 8 \
   --keep-tmp-files

A detailed description of all command options can be found on the :ref:`Usage <pvacvector_run>` page.
