.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

Output Files
============

pVACvector will create the following output files:

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
   * - ``vector.png``
     - A PNG visualization of the above result.
   * - ``<sample_name>_results.dna.fa``
     - The final output file with the backtranslated DNA sequences of the included peptides and best spacers in the optimal order.
   * - ``without_<peptide_id(s)>`` (directory
     - If no solution was found with spacers and clipping peptides, these
       directories contain result files for partial solutions obtained by removing
       peptide(s) as indicated in each directory name. These directories
       in turn contain a ``junctions.tsv`` file of all of the junctions
       remaining after the respective peptide(s) were removed, as well as a
       ``<sample_name>_results.fa``, a ``<sample_name>_results.dna.fa``, and a
       ``vector.png`` file if a solution was found by removing the indicated
       peptide(s).

.. figure:: ../images/vector.jpg
   :align: center
   :alt: pVACvector result visualization example

   pVACvector result visualization example
