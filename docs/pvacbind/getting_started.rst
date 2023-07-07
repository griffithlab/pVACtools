.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

Getting Started
---------------

pVACbind provides a set of example data to show the expected format of input and output files.
You can download the data set by running the ``pvacbind download_example_data`` :ref:`command <pvacbind_example_data>`.

The example data output can be reproduced by running the following command:

.. code-block:: none

   pvacbind run \
   <example_data_dir>/input.fasta \
   Test \
   HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
   all \
   <output_dir> \
   -e1 8,9,10 \
   -e2 15

A detailed description of all command options can be found on the :ref:`Usage <pvacbind_run>` page.
