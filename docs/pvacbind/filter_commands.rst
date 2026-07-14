.. image:: ../images/pVACbind_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACbind logo

.. _pvacbind_filter_commands:

Filtering Commands
=============================

pVACbind currently offers two filters: a binding filter and a top score filter.

These filters are always run automatically as part
of the pVACbind pipeline using default cutoffs.

All filters can also be run manually on the filtered.tsv file to narrow the results down further,
or they can be run on the all_epitopes.tsv file to apply different filtering thresholds.

The binding filter is used to remove neoantigen candidates that do not meet desired peptide:MHC binding criteria.
The top score filter is used to select the most promising peptide candidate for each variant. 
Multiple candidate peptides from a single somatic variant can be caused by multiple peptide lengths, registers, HLA alleles,
and transcript annotations.

Further details on each of these filters is provided below.

.. note::

   The default values for filtering thresholds are suggestions only. While they are based on review of the literature
   and consultation with our clinical and immunology colleagues, your specific use case will determine the appropriate values.

Binding Filter
--------------

.. program-output:: pvacbind binding_filter -h

The binding filter filters out variants that don't pass the chosen binding threshold.
The user can chose whether to apply this filter to the ``lowest`` or the ``median`` binding
affinity score by setting the ``--top-score-metric`` flag. The ``lowest`` binding
affinity score is recorded in the ``Best IC50 Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median IC50 Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.
An additional ``--top-score-metric2`` flag allows the user to choose whether to use IC50 or
Percentile scores. By default, IC50 is used.

When the ``--allele-specific-binding-thresholds`` flag is set, binding cutoffs specific to each
prediction's HLA allele are used instead of the value set via the ``--binding-threshold`` parameters.
For HLA alleles where no allele-specific binding threshold is available, the
binding threshold is used as a fallback. Alleles with allele-specific
threshold as well as the value of those thresholds can be printed by executing
the ``pvacbind allele_specific_cutoffs`` command.

In addition to being able to filter on the IC50 score columns, the binding
filter also offers the ability to filter on the percentile score using the
``--percentile-threshold`` parameter. When the ``--top-score-metric`` is set
to ``lowest``, this threshold is applied to the ``Best Percentile`` column. When
it is set to ``median``, the threshold is applied to the ``Median
Percentile`` column.

When the ``--percentile-threshold`` flag is set, the candidate inclusion strategy can be
specified by using the ``--percentile-threshold-strategy`` parameter. The parameter has two
options ``conservative`` (default) and ``exploratory``. The 'conservative' option requires a candidate 
to pass BOTH the binding threshold and percentile threshold, while the 'exploratory' option requires
a candidate to pass EITHER the binding threshold or percentile threshold.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Top Score Filter
----------------

.. program-output:: pvacbind top_score_filter -h

This filter picks the top epitope for each peptide in the original input FASTA.
Epitopes with the same Index are identified as coming from the same peptide.

For each peptide the best epitope is determined as follows:

- Pick the entries with no Problematic Positions.
- For the remaining entries, calculate a rank for all the metrics specified
  via the ``--top-score-metric2`` parameter and sum them. Whether the lowest or median value
  is considered for each metric is controlled by the ``--top-score-metric`` parameter.
  Sort the remaining entries on this sum rank followed by the rank of the first
  ``--top-score-metric2`` specified (to break
  any ties in the sum rank). Select the highest sorted entry.

The selected top epitopes for each peptide are then sorted as follows:

.. list-table::
   :header-rows: 1

   * - Sort Criteria
     - Sort Order
   * - Sum of the ascending ranks of
       the metrics selected via the ``--top-score-metric2`` parameter (possible values:
       ``IC50 MT``, ``%ile MT``, ``IC50 %ile MT``, ``Pres %ile MT``; default: ``IC50 MT``,
       ``%ile MT``).
     - Ascending sum rank
   * - First metric specified in the ``--top-score-metric2`` as a tie breaker
       for identical sum ranks
     - Ascending rank
   * - ``Index`` column
     - Alphabetical

Aggregate Report Filter
-----------------------

.. program-output:: pvacbind aggregate_report_filter -h

This command filters the aggregate report to only those variants matching the
specified ``--include-tiers`` (default:Pass).
