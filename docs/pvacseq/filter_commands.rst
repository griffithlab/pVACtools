.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Filtering Commands
=============================

pVACseq currently offers three filters: a binding filter, a coverage filter,
and a top score filter.

These filters are always run automatically as part
of the pVACseq pipeline using default cutoffs.

All filters can also be run manually on the filtered.tsv file to narrow the results down further
or they can be run on the all_epitopes.tsv file to apply different filtering
thresholds.

Binding Filter
--------------

.. program-output:: pvacseq binding_filter -h

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

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Coverage Filter
---------------

.. program-output:: pvacseq coverage_filter -h

.. .. argparse::
    :module: lib.coverage_filter
    :func: define_parser
    :prog: pvacseq coverage_filter

If the input VCF contains readcount and/or expression annotations,
then the coverage filter
can be run again on the filtered.tsv report file to narrow down the results even further.
You can also run this filter again on the all_eptiopes.tsv report file to
apply different cutoffs.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Top Score Filter
----------------

.. program-output:: pvacseq top_score_filter -h

This filter picks the top epitope for a variant. By default the
``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median MT Score`` column and pick the epitope with the lowest
median mutant ic50 score for each variant. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best MT Score`` column is instead used to
make this determination.
