.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _filter_commands:

Filtering Commands
=============================

pVACseq currently offers four filters: a binding filter, a coverage filter,
a transcript filter, and a top score filter.

These filters are always run automatically as part
of the pVACseq pipeline using default cutoffs.

All filters can also be run manually on the filtered.tsv file to narrow the results down further,
or they can be run on the all_epitopes.tsv file to apply different filtering thresholds.

The binding filter is used to remove neoantigen candidates that do not meet desired peptide:MHC binding criteria.
The coverage filter is used to remove variants that do not meet desired read count and VAF criteria (in normal DNA
and tumor DNA/RNA). The transcript filter is used to remove variant annotations based on low quality
transcript annotations. The top score filter is used to select the most promising peptide candidate for each variant.
Multiple candidate peptides from a single somatic variant can be caused by multiple peptide lengths, registers, HLA alleles,
and transcript annotations.

Further details on each of these filters is provided below.

.. note::

   The default values for filtering thresholds are suggestions only. While they are based on review of the literature
   and consultation with our clinical and immunology colleagues, your specific use case will determine the appropriate values.

Binding Filter
--------------

.. program-output:: pvacseq binding_filter -h

.. .. argparse::
    :module: lib.binding_filter
    :func: define_parser
    :prog: pvacseq binding_filter

The binding filter removes variants that don't pass the chosen binding threshold.
The user can chose whether to apply this filter to the ``lowest`` or the ``median`` binding
affinity score by setting the ``--top-score-metric`` flag. The ``lowest`` binding
affinity score is recorded in the ``Best MT IC50 Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median MT IC50 Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.
An additional ``--top-score-metric2`` flag allows the user to choose whether to use IC50 or
Percentile scores. By default, IC50 is used.

When the ``--allele-specific-binding-thresholds`` flag is set, binding cutoffs specific to each
prediction's HLA allele are used instead of the value set via the ``--binding-threshold`` parameters.
For HLA alleles where no allele-specific binding threshold is available, the
binding threshold is used as a fallback. Alleles with allele-specific
threshold as well as the value of those thresholds can be printed by executing
the ``pvacseq allele_specific_cutoffs`` command.

The binding filter also offers the option to filter on ``Fold Change`` columns, which contain
the ratio of the MT score to the WT Score. This option can be activated by setting the
``--minimum-fold-change`` threshold (to require that the mutant peptide is a better binder 
than the corresponding wild type peptide). If the ``--top-score-metric`` option is set to ``lowest``, 
the ``Corresponding Fold Change`` column will be used (``Corresponding WT IC50 Score``/``Best MT IC50 Score``).
If the ``--top-score-metric`` option is set to ``median``, the ``Median Fold Change`` column
will be used (``Median WT IC50 Score``/``Median MT IC50 Score``).

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

.. program-output:: pvacseq coverage_filter -h

.. .. argparse::
    :module: lib.coverage_filter
    :func: define_parser
    :prog: pvacseq coverage_filter

If the input VCF contains readcount and/or expression annotations, then the coverage filter
can be run again on the filtered.tsv report file to narrow down the results even further.
You can also run this filter again on the all_epitopes.tsv report file to apply different cutoffs. 

The general goals of these filters are to limit variants for neoepitope prediction to those 
with good read support and/or remove possible sub-clonal variants. In some cases the input 
VCF may have already been filtered in this fashion. This filter also allows for removal of
variants that do not have sufficient evidence of RNA expression.

For more details on how to prepare input VCFs that contain all of these annotations, refer to 
the :ref:`prerequisites_label` section for more information.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Transcript Filter
-----------------

.. program-output:: pvacseq transcript_filter -h

This filter is used to eliminate variant annotations based on poorly-supported transcripts. This assessed
based on whether the transcript is the MANE Select transcripts, whether it is
the canonical transcript or whether the transcript support level (TSL) meets the
``--maximum-transcript-support-level`` cutoff. The
``--transcript-prioritizatio-strategy`` parameter controls which ones of these three
criteria are considered. A neoantigen candidate passes this filter if its
transcript passes at least one of the specified criteria.

Transcript with a TSL of ``Not Supported`` will pass the TSL criteria. These values occur if VEP was run
without the ``--tsl`` flag or if data is aligned to GRCh37 or older.

Top Score Filter
----------------

.. program-output:: pvacseq top_score_filter -h

This filter picks the top epitope for a variant. Epitopes with the same
Chromosome - Start - Stop - Reference - Variant are identified as coming from
the same variant.

In order to account for different splice sites among the transcripts of a
variant that would lead to different peptides, this filter also takes into
account the different transcripts returned by VEP and bins the ones resulting
in the same set of epitopes together into a transcript set. For each transcript
set the filter will return the top epitope similar to how the Best Peptide is
determined in the :ref:`aggregated report <aggregated>`:

- If ``--allow-inclomplete-transcripts`` flag is set, pick the entries without
  a Transcript CDS Flags set.
- Of the remaining entries, pick the entries where the Biotype is ``protein_coding``.
- Of the remaining entries, pick the entries that pass at least one of the transcript criteria selected in the
  ``--transcript-prioritization-strategy`` taking into consideration the
  ``--maximum-transcript-support-level`` if ``tsl`` is one of the selected
  criteria.
- Of the remaining entries, pick the entries with no Problematic Positions.
- Of the remaining entries, pick the ones passing the Anchor Criteria (see
  details below)
- For the remaining entries, calculate a rank for all the metrics specified
  via the ``--top-score-metric2`` parameter and sum them. Whether the lowest or median value
  is considered for each metric is controlled by the ``--top-score-metric`` parameter.
  Sort the remaining entries on this sum rank followed by the rank of the first
  ``--top-score-metric2`` specified (to break
  any ties in the sum rank), MANE Select status, Canonical status,
  Transcript Support Level, Transcript Length, and Transcript
  Expression. Select the highest sorted entry.

The selected top epitopes for each transcript set are then sorted as follows:

.. list-table::
   :header-rows: 1

   * - Sort Criteria
     - Sort Order
   * - Sum of ascending ranks of ``Gene Expression`` and the ascending ranks of
       the metrics selected via the ``--top-score-metric2`` parameter (possible values:
       ``IC50 MT``, ``%ile MT``, ``IC50 %ile MT``, ``Pres %ile MT``; default: ``IC50 MT``,
       ``%ile MT``).
     - Ascending sum rank
   * - First metric specified in the ``--top-score-metric2`` as a tie breaker
       for identical sum ranks
     - Ascending rank
   * - ``Gene Name`` column
     - Alphabetical
   * - ``Mutation`` column
     - Alphabetical

**Anchor Criteria**

This criteria is failed if all mutated amino acids of the entry (``Mutation
Position``) are at an anchor position and the WT peptide has good binding
``(Best/Median WT IC50 Score < binding_threshold)``.

When the ``--allele-specific-binding-thresholds`` flag is set, binding cutoffs specific to each
prediction's HLA allele are used instead of the value set via the ``--binding-threshold`` parameters.
For HLA alleles where no allele-specific binding threshold is available, the
binding threshold is used as a fallback. Alleles with allele-specific
threshold as well as the value of those thresholds can be printed by executing
the ``pvacseq allele_specific_cutoffs`` command.

**Additional Considerations**

It is important to note that there are several reasons why a particular variant can lead to multiple peptides
with different predicted binding affinities. The following can result in multiple peptides and/or binding predictions for a single
variant:

#. Different epitope lengths: specifying multiple epitope lengths results in similar but non-identical epitope sequences for each variant (e.g. KLPEPCPS, KLPEPCPST, KLPEPCPSTT, KLPEPCPSTTP).
#. Different registers: pVACseq will test epitopes where the mutation is in every position (e.g. EPCPSTTP, PEPCPSTT, LPEPCPST, KLPEPCPS, ...).
#. Different transcripts: in some case the peptide sequence surrounding a variant will depend on the reference transcript sequence, particularly if there are alternative splice sites near the variant position.
#. Different HLA alleles: the HLA allele that produces the best predicted binding affinity is chosen.
#. A homozygous somatic variant with heterozygous proximal variants nearby may produce multiple different peptides.

The significance of choosing a single representative peptide can depend on your experimental or clinical aims.
For example, if you are planning to use short peptide sequences exactly as they were assessed 
for binding affinity in pVACseq (e.g. specific 9-mers for in vitro experimental validation or perhaps a dendritic cell vaccine delivery 
approach) then the selection of a specific peptide from the possibilities caused by different lengths, registers, etc. 
is very important. In some cases you may wish to consider more criteria beyond which of these candidates has the best 
predicted binding affinity and gets chosen by the Top Score Filter. 

On the other hand, if you plan to use synthetic long peptides (SLPs) or encode your candidates in a DNA vector, you will likely include 
flanking amino acids. This means that you often get a lot of the different short peptides that correspond to slightly different lengths or 
registers within the longer containing sequence. In this scenario, pVACseq's choice of a single candidate peptide by the Top Score Filter 
isn't actually that critical in the sense of losing other good candidates, because you may get them all anyway.

One important exception to this is the rare case where the same variant leads to different peptides in different transcripts (due to different splice site usage).
If multiple transcripts are expressed and 
lead to distinct peptides, you may want to include both in your final list of candidates.
The top score filter supports this case, as described above.
This assumes you did not start with only a single transcript
model for each gene (e.g. using the ``--pick`` option in VEP) and also that if you are requiring transcripts with TSL=1 that there
are multiple qualifying transcripts that lead to different peptide sequences at the site of the variant. This will be fairly rare.
Even though most genes have alternative transcripts, they often have only subtle differences in open reading frame and overall
protein sequence, and only differences within the window that would influence a neoantigen candidate are consequential here.

Aggregate Report Filter
-----------------------

.. program-output:: pvacseq aggregate_report_filter -h

This command filters the aggregate report and its corresponding metrics.json
file to only those variants matching the specified ``--include-tiers`` (default:
Pass).

This filter may be useful for high tumor mutation burden cases when the
metrics.json file size exceeds the pVACview upload file size limit or
slows down the pVACview application.
