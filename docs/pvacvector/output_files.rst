.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Output Files
============

The pVACseq pipeline will write its results in separate folders depending on
which prediction algorithms were chosen:

- ``MHC_Class_I``: for MHC class I prediction algorithms
- ``MHC_Class_II``: for MHC class II prediction algorithms
- ``combined``: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoepitope predictions from both

Each folder will contain the same list of output files (listed in the order
created):

.. list-table::
   :header-rows: 1

   * - File Name
     - Description
   * - ``vector_input.fa`` (optional)
     - An intermediate file with vaccine peptide sequences created from the epitopes in a pVACseq output file.
   * - 0...n (directory)
     - One numbered directory for each iteration of clipping necessary. Each
       one contains subdirectories for the spacers tested which in turn contain
       prediction information and a ``junctions.tsv`` intermediate result file for the attempt iteration.
   * - ``<sample_name>_results.fa``
     - The final output file with the peptide sequences and best spacers in the optimal order.
   * - ``junctions.tsv``
     - A tab-separated file listing all of the valid junctions found by pVACvector including spacer and clipping information.
   * - ``vector.jpg``
     - A JPEG visualization of the above result.
   * - ``<sample_name>_results.dna.fa``
     - The final output file with the backtranslated DNA sequences of the included peptides and best spacers in the optimal order.

.. figure:: ../images/vector.jpg
   :align: center
   :alt: pVACvector result visualization example

   pVACvector result visualization example
