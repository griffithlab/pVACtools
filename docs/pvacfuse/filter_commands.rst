.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Filtering Commands
=============================

pVACfuse currently offers a binding filter.

The binding filter is always run automatically as part of the pVACfuse pipeline. It can also be run manually to narrow the final results down further.

Binding Filter
--------------

.. topic:: For usage instructions run

   ``pvacfuse binding_filter --help``

.. .. argparse::
    :module: lib.binding_filter
    :func: define_parser
    :prog: pvacfuse binding_filter

The binding filter filters out variants that don't pass the chosen binding threshold. The user can chose whether to apply this filter to the "lowest" or the "median" binding affinity score. The "lowest" binding affinity score is recorded in the "Best MT Score" column and represents the lowest ic50 score of all prediction algorithms that were picked during the previous pVACfuse run. The "median" binding affinity score is recorded in the "Median MT Score" column and corresponds to the median ic50 score of all prediction algorithms used to create the report.

The binding filter also offers the option to filter on Fold Change columns, which contain the ratio of the MT score to the WT Score. If the binding filter is set to "best", the "Corresponding Fold Change" column will be used. ("Corresponding WT Score"/"Best MT Score"). If the binding filter is set to "median", the "Median Fold Change" column will be used ("Median WT Score"/"Median MT Score").

.. Coverage Filter
 ---------------

.. .. topic:: For usage instructions run  
  .. ``pvacfuse coverage_filter --help``

.. .. argparse::
    :module: lib.coverage_filter
    :func: define_parser
    :prog: pvacseq coverage_filter

.. If a pVACfuse process has been run with bam-readcount or Cufflinks input files then the coverage_filter can be run again on the final report file to narrow down the results even further.

.. If no additional coverage input files have been provided to the main pVACfuse run then this information would need to be manually added to the report in order to run this filter.

Top Score Filter
----------------

.. topic:: For usage instructions run

   ``pvacfuse top_score_filter --help``

This filter picks the top epitope for a variant. By default the
``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median MT Score`` column and pick the epitope with the lowest
median mutant ic50 score for each variant. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best MT Score`` column is instead used to
make this determination.
