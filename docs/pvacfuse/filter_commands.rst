.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

.. _pvacfuse_filter_commands:

Filtering Commands
=============================

pVACfuse currently offers three filters: a binding filter, a coverage filter,
and a top score filter.

All filters are run automatically as part of the pVACfuse pipeline.

All filters can also be run manually to narrow the final results down further 
or to redefine the filters entirely and produce a new candidate list from the 
all_epitopes.tsv file.

.. note::

   The default values for filtering thresholds are suggestions only. While they are based on review of the literature
   and consultation with our clinical and immunology colleagues, your specific use case will determine the appropriate values.

Binding Filter
--------------

.. program-output:: pvacfuse binding_filter -h

.. .. argparse::
    :module: lib.binding_filter
    :func: define_parser
    :prog: pvacfuse binding_filter

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
the ``pvacfuse allele_specific_cutoffs`` command.

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

Coverage Filter
---------------

.. program-output:: pvacfuse coverage_filter -h

If a pVACfuse process has been run with Arriba data, Read Support information will be available.
If AGFusion data was used an input, a STAR-Fusion file will have needed to be
provided in the run in order to make Read Support and Expression information available.

The coverage filter
can be run again on the filtered.tsv report file to narrow down the results even further.
You can also run this filter on the all_epitopes.tsv report file to apply different cutoffs.

The general goals of this filter is to limit variants for neoepitope prediction to those 
with good read support. In some cases the input data may have already been filtered in this fashion.
This filter also allows for removal of variants that do not have sufficient evidence of RNA expression.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Top Score Filter
----------------

.. program-output:: pvacfuse top_score_filter -h

This filter picks the top epitope for a fusion. Epitopes with the same
Chromosome - Start - Stop are identified as coming from the same fusion.

In order to account for different splice sites among the transcripts of a
fusion that would lead to different peptides, this filter also takes into
account the different transcripts returned by AGFusion/Arriba and bins the
ones resulting in the same set of epitopes together into a transcript set.
For each transcript set the filter will return the top epitope as follows:

- Pick the entries with no Problematic Positions.
- For the remaining entries, calculate a rank for all the metrics specified
  via the ``--top-score-metric2`` parameter and sum them. Whether the lowest
  or median value is considered for each metric is controlled by the
  ``--top-score-metric`` parameter. Sort the remaining entries on this sum
  rank followed by the rank of the first ``top-score-metric2``  specified
  (to break any ties in the sum rank), and Expression. Select the highest
  sorted entry.

The selected top epitopes for each transcript set are then sorted as follows:

.. list-table::
   :header-rows: 1

   * - Sort Criteria
     - Sort Order
   * - Sum of ascending ranks of ``Exprission`` and the ascending ranks of
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

.. program-output:: pvacfuse aggregate_report_filter -h

This command filters the aggregate report to only those variants matching the
specified ``--include-tiers`` (default:Pass).
