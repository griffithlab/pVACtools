.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Filtering Commands
=============================

pVACfuse currently offers two filters: a binding filter
and a top score filter.

The binding filter and top score filter are always run automatically as part
of the pVACfuse pipeline.

All filters can also be run manually to narrow the final results down further 
or to redefine the filters entirely and produce a new candidate list from the 
all_epitopes.tsv file.

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
affinity score is recorded in the ``Best MT Score`` column and represents the lowest
ic50 score of all prediction algorithms that were picked during the previous pVACseq run.
The ``median`` binding affinity score is recorded in the ``Median MT Score`` column and
corresponds to the median ic50 score of all prediction algorithms used to create the report.
Be default, the binding filter runs on the ``median`` binding affinity.

By default, entries with ``NA`` values will be included in the output. This
behavior can be turned off by using the ``--exclude-NAs`` flag.

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

.. program-output:: pvacfuse top_score_filter -h

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
