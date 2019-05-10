.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _filter_commands:

Filtering Commands
=============================

pVACseq currently offers four filters: a binding filter, a coverage filter,
a transcript support level filter, and a top score filter.

These filters are always run automatically as part
of the pVACseq pipeline using default cutoffs.

All filters can also be run manually on the filtered.tsv file to narrow the results down further,
or they can be run on the all_epitopes.tsv file to apply different filtering thresholds.

The binding filter is used to remove neoantigen candidates that do not meet desired peptide:MHC binding criteria.
The coverage filter is used to remove variants that do not meet desired read count and VAF criteria (in normal DNA
and tumor DNA/RNA). The transcript support level filter is used to remove variant annotations based on low quality
transcript annotations. The top score filter is used to select the most promising peptide candidate for each variant. 
Multiple candidate peptides from a single somatic variant can be caused by multiple peptide lengths, registers, HLA alleles,
and transcript annotations.

Further details on each of these filters is provided below.

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
affinity score is recorded in the ``Best MT Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median MT Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.

The binding filter also offers the option to filter on ``Fold Change`` columns, which contain
the ratio of the MT score to the WT Score. This option can be activated by setting the
``--minimum-fold-change`` threshold (to require that the mutant peptide is a better binder 
than the corresponding wild type peptide). If the ``--top-score-metric`` option is set to ``lowest``, 
the ``Corresponding Fold Change`` column will be used (``Corresponding WT Score``/``Best MT Score``).
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

Transcript Support Level Filter
-------------------------------

.. program-output:: pvacseq transcript_support_level_filter -h

This filter is used to eliminate variant annotations based on poorly-supported transcripts. By default,
only transcripts with a `transcript support level (TSL) <https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl>`_
of <=1 are kept. This threshold can be adjusted using the ``--maximum-transcript-support-level``
parameter. 

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

Top Score Filter
----------------

.. program-output:: pvacseq top_score_filter -h

This filter picks the top epitope for a variant. Epitopes with the same
Chromosome - Start - Stop - Reference - Variant are identified as coming from
the same variant.

By default the
``--top-score-metric`` option is set to ``median`` which will apply this
filter to the ``Median MT Score`` column and pick the epitope with the lowest
median mutant ic50 score for each variant. If the ``--top-score-metric``
option is set to ``lowest``, the ``Best MT Score`` column is instead used to
make this determination.

If there are multiple top epitopes for a variant with the same ic50 score, the
first one is chosen.

It is important to note that there are several reasons why a particular variant can lead to multiple peptides
with different predicted binding affinities. The following can result in multiple peptides and/or binding predictions for a single
variant:

1. Different peptide lengths. e.g. if you specify multiple lengths, you get similar but non-idential peptide sequences for each 
variant (e.g. KLPEPCPS, KLPEPCPST, KLPEPCPSTT, KLPEPCPSTTP).
2. Different registers. pvacseq will test peptides where the mutation is in every postion (e.g. EPCPSTTP, PEPCPSTT, LPEPCPST, KLPEPCPS, ...).
3. Different transcripts. In some case the peptide sequence surrounding a variant will depend on the reference transcript sequence, particularly
if there are alternative splice sites near the variant position.
4. Different HLA alleles. The HLA allele that produces the best predicted binding affinity is chosen.
5. A homozygous somatic variant with heterozygous proximal variants nearby may produce multiple different peptides.

The significance of choosing a single representative peptide can depend on your experimental or clinical aims.
For example, if you are planning to use short peptide sequences exactly as they were assessed 
for binding affinity in pvacseq (e.g. specific 9-mers for in vitro experimental validation or perhaps a dendritic cell vaccine delivery 
approach) then the selection of a specific peptide from the possibilities caused by different lengths, registers, etc. 
is very important. In some cases you may wish to consider more criteria beyond which of these candidates has the best 
predicted binding affinity and gets chosen by the Top Score Filter. 

On the other hand, if you plan to use synthetic long peptides (SLPs) or encode your candidates in a DNA vector, you will likely include 
flanking amino acids. This means that you often get a lot of the different short peptides that correspond to slightly different lengths or 
registers within the longer containing sequence. In this scenario, pvacseq's choice of a single candidate peptide by the Top Score Filter 
isn't actually that critical in the sense of losing other good candidates, because you may get them all anyway.

One important exception to this is the rare case where the same variant leads to different peptides in different transcripts 
(due to different splice site usage). In this case, pvacseq still selects the top candidate as described above but including flanking 
amino acid sequence won't help because you are only doing this from a single transcript. If multiple transcripts are expressed and 
lead to distinct peptides, you may want to include both in your final list of candidates. The Top Score Filter doesn't yet support 
this but the unfiltered results will still contain these candidates. This assumes you did not start with only a single transcript
model for each gene (e.g. using the --pick option in ensembl) and also that if you are requiring transcripts with TSL=1 that there
are multiple qualifying transcripts that lead to different peptide sequences at the site of the variant. This will be fairly rare.
Even though most genes have alternative transcripts, they often have only subtle differences in open reading frame and overall
protein sequence, and only differences within the window that would influence a neoantigen candidate are consequential here.

