.. image:: ../images/pVACsplice_logo_trans-bg_v4b.png
    :align: right
    :alt: pVACsplice logo
    :width: 175px

.. _pvacsplice_filter_commands:

Filtering Commands
==================

pVACsplice currently offers four filters: a binding filter, a coverage filter,
a transcript support level filter, and a top score filter.

These filters are always run automatically as part
of the pVACsplice pipeline using default cutoffs.

All filters can also be run manually on the filtered.tsv file to narrow the results down further,
or they can be run on the all_epitopes.tsv file to apply different filtering thresholds.

The binding filter is used to remove neoantigen candidates that do not meet desired peptide:MHC binding criteria.
The coverage filter is used to remove variants that do not meet desired read count and VAF criteria (in normal DNA
and tumor DNA/RNA). The transcript support level filter is used to remove variant annotations based on low quality
transcript annotations. The top score filter is used to select the most promising peptide candidate for each variant.
Multiple candidate peptides from a single somatic variant can be caused by multiple peptide lengths, registers, HLA alleles,
and transcript annotations.

Further details on each of these filters is provided below.

.. note::

   The default values for filtering thresholds are suggestions only. While they are based on review of the literature
   and consultation with our clinical and immunology colleagues, your specific use case will determine the appropriate values.

Binding Filter
--------------

.. program-output:: pvacsplice binding_filter -h

The binding filter removes variants that don't pass the chosen binding threshold.
The user can chose whether to apply this filter to the ``lowest`` or the ``median`` binding
affinity score by setting the ``--top-score-metric`` flag. The ``lowest`` binding
affinity score is recorded in the ``Best MT IC50 Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median MT IC50 Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.

When the ``--allele-specific-binding-thresholds`` flag is set, binding cutoffs specific to each
prediction's HLA allele are used instead of the value set via the ``--binding-threshold`` parameters.
For HLA alleles where no allele-specific binding threshold is available, the
binding threshold is used as a fallback. Alleles with allele-specific
threshold as well as the value of those thresholds can be printed by executing
the ``pvacsplice allele_specific_cutoffs`` command.

In addition to being able to filter on the IC50 score columns, the binding
filter also offers the ability to filter on the percentile score using the
``--percentile-threshold`` parameter. When the ``--top-score-metric`` is set
to ``lowest``, this threshold is applied to the ``Best MT Percentile`` column. When
it is set to ``median``, the threshold is applied to the ``Median MT
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

.. program-output:: pvacsplice coverage_filter -h

If the pVACsplice input VCF contains readcount and/or expression annotations, then the coverage filter
can be run again on the filtered.tsv report file to narrow down the results even further.
You can also run this filter again on the all_epitopes.tsv report file to apply different cutoffs.

The general goals of these filters are to limit variants for neoepitope prediction to those
with good read support and/or remove possible sub-clonal variants. In some cases the input
VCF may have already been filtered in this fashion. This filter also allows for removal of
variants that do not have sufficient evidence of RNA expression.

For more details on how to prepare input VCFs that contain all of these annotations, refer to
the :ref:`pvacsplice_prerequisites_label` section for more information.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Transcript Support Level Filter
-------------------------------

.. program-output:: pvacsplice transcript_support_level_filter -h

This filter is used to eliminate variant annotations based on poorly-supported transcripts. By default,
only transcripts with a `transcript support level (TSL) <https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl>`_
of <=1 are kept. This threshold can be adjusted using the ``--maximum-transcript-support-level``
parameter.

By default, entries with ``Not Supported`` values will be included in the output. These occur if VEP was run
without the ``--tsl`` flag or if data is aligned to GRCh37 or older.

Top Score Filter
----------------

.. program-output:: pvacsplice top_score_filter -h

This filter picks the top epitope for each splice site variant. The top epitope is
determined by first selecting epitopes with no Problematic Positions
and among those selecting the one with lowest median/best MT IC50 score for
each splice site variant

By default the ``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median MT IC50 Score`` column. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best MT IC50  Score`` column is used
instead.
