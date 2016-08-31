Filtering Commands
=============================

pVAC-Seq currently offers two filters: a binding filter and a coverage filter.

The binding filter is always run automatically as part of the pVAC-Seq pipeline but can be run again on the final output file if further filtering is desired.

The coverage filter is run automatically if bam-readcount or cufflinks file are provided as additional input files to a pVAC-Seq run but can also be run again on the final report file to narrow the results down further.

Binding Filter
--------------

.. argparse::
    :module: lib.binding_filter
    :func: define_parser
    :prog: pvacseq binding_filter

    The binding filter filters out variants that don't pass the chosen binding threshold. The user can chose whether to apply this filter to the "lowest" or the "median" binding affinity score. The "lowest" binding affinity score is recorded in the "Best MT Score" column and represents the lowest ic50 score between all prediction algorithms that were picked during the previous pVAC-Seq run. The "median" binding affinity score is recorded in the "Median MT Score" column and corresponds to the median ic50 score between all prediction algorithms used to create the report.

    The binding filter also offers the option to filter on the "Corresponding Fold Change" column, which is the ratio of the "Best MT Score" to the "Corresponsing WT Score".

Coverage Filter
---------------

.. argparse::
    :module: lib.coverage_filter
    :func: define_parser
    :prog: pvacseq coverage_filter

    If a pVAC-Seq process has been run with bam-readcount or Cufflinks input files then the coverage_filter can be run again on the final report file to narrow down the results even further.

    If no additional coverage input files have been provided to the main pVAC-Seq run then this information would need to be manually added to the report in order to run this filter.
