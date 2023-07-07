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

This filter picks the top epitope for a variant. Epitopes with the same
Chromosome - Start - Stop - Reference - Variant are identified as coming from
the same variant.

In order to account for different splice sites among the transcripts of a
variant that would lead to different peptides, this filter also takes into
account the different transcripts returned by AGFusion/Arriba and will return
the top epitope for each transcript if they are non-identical. If the
resulting list of top epitopes for the transcripts of a variant is identical,
the epitope for the transcript with the highest expression is returned. If
this information is not available, the transcript with the lowest Ensembl ID is returned.

By default the
``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median IC50 Score`` column and pick the epitope with the lowest
median mutant ic50 score for each variant. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best IC50 Score`` column is instead used to
make this determination.
