.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Filtering Commands
=============================

pVACseq currently offers three filters: a binding filter, a coverage filter,
and a top score filter.

The binding filter is always run automatically as part of the pVACseq pipeline.
The coverage filter is run automatically if bam-readcount or cufflinks file are
proAvided as additional input files to a pVACseq run. The top score filter is
run if the ``--top-result-per-mutation`` flag is set.

All filters can also be run manually to narrow the final results down further.

Binding Filter
--------------

.. topic:: For usage instructions run

   ``pvacseq binding_filter --help``

.. .. argparse::
    :module: lib.binding_filter
    :func: define_parser
    :prog: pvacseq binding_filter

The binding filter filters out variants that don't pass the chosen binding threshold.
The user can chose whether to apply this filter to the ``lowest`` or the ``median`` binding
affinity score by setting the ``--top-score-metric`` flag. The ``lowest`` binding
affinity score is recorded in the ``Best MT Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median MT Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.

The binding filter also offers the option to filter on ``Fold Change`` columns, which contain
the ratio of the MT score to the WT Score. This option can be activated by setting the
``--minimum-fold-change`` threshold. If the ``--top-score-metric`` option is set to ``lowest``, the
``Corresponding Fold Change`` column will be used (``Corresponding WT Score``/``Best MT Score``).
If the ``--top-score-metric`` option is set to ``median``, the ``Median Fold Change`` column
will be used (``Median WT Score``/``Median MT Score``).

Coverage Filter
---------------

.. topic:: For usage instructions run

   ``pvacseq coverage_filter --help``

.. .. argparse::
    :module: lib.coverage_filter
    :func: define_parser
    :prog: pvacseq coverage_filter

If a pVACseq process has been run with bam-readcount or Cufflinks input files then the coverage filter
can be run again on the final report file to narrow down the results even further.

If no additional coverage input files have been provided to the main pVACseq run then this information
would need to be manually added to the report in order to run this filter
using the appropriate headers. Columns available for this filter are ``Tumor DNA Depth``, ``Tumor DNA VAF``,
``Tumor RNA Depth``, ``Tumor RNA VAF``, ``Normal Depth``, ``Normal VAF``, ``Gene Expression``, ``Transcript Expression``.

Top Score Filter
----------------

.. topic:: For usage instructions run

   ``pvacseq top_score_filter --help``

This filter picks the top epitope for a variant. By default the
``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median MT Score`` column and pick the epitope with the lowest
median mutant ic50 score for each variant. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best MT Score`` column is instead used to
make this determination.
