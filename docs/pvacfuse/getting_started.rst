.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Getting Started
---------------

pVACfuse provides a set of example data to show the expected format of input and output files. 
You can download the data set by running the ``pvacfuse download_example_data`` :ref:`command <pvacfuse_example_data>`.

The AGFusion example data output can be reproduced by running the following command:

.. code-block:: none

   pvacfuse run \
   <example_data_dir>/agfusion/ \
   Test \
   HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
   all \
   <output_dir> \
   -e1 8,9,10 \
   -e2 15

The Arriba example data output can be reproduced by running the following
command:

.. code-block:: none

   pvacfuse run \
   <example_data_dir>/arriba_fusions.tsv \
   Test \
   HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
   all \
   <output_dir> \
   -e1 8,9,10 \
   -e2 15

A detailed description of all command options can be found on the :ref:`Usage <pvacfuse_run>` page.
