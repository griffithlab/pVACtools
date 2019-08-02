.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Getting Started
---------------

pVACfuse provides a set of example data to show the expected format of input and output files. 
You can download the data set by running the ``pvacfuse download_example_data`` :ref:`command <pvacfuse_example_data>`.

There are two option as to how to run pVACfuse. It accepts either a
INTEGRATE-neo output bedpe file or a AGFusion output directory.

The following command is an example for how to run pVACfuse with an
INTEGRATE-neo bedpe file and will regenerate the
``results_from_integrate_neo`` example data:

.. code-block:: none

   pvacfuse run \
   <example_data_dir>/fusions.bedpe.annot \
   Test \
   HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
   MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
   <output_dir> \
   -e 8,9,10

The ``results_from_agfusion`` example data can be regenerated like so:

.. code-block:: none

   pvacfuse run \
   <example_data_dir>/agfusion/ \
   Test \
   HLA-A*02:01,HLA-B*35:01,DRB1*11:01 \
   MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
   <output_dir> \
   -e 8,9,10

A detailed description of all command options can be found on the :ref:`Usage <pvacfuse_run>` page.
