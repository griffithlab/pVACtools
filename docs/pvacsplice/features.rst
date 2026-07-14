.. image:: ../images/pVACsplice_logo_trans-bg_v4b.png
    :align: right
    :alt: pVACsplice logo
    :width: 175px

Features
========

**Splice Site Analysis**

pVACsplice offers epitope binding predictions for splice site variants
predicted by RegTools.

**No local install of epitope prediction software needed**

pVACsplice utilizes the IEDB RESTful web interface. This means that none of the underlying prediction software, like NetMHC, needs to be installed locally.

.. warning::
   We only recommend using the RESTful API for small requests. If you use the
   RESTful API to process large VCFs or to make predictions for many alleles,
   epitope lengths, or prediction algorithms, you might overload their system.
   This can result in the blacklisting of your IP address by IEDB, causing
   403 errors when trying to use the RESTful API. In that case please open
   a ticket with `IEDB support <http://help.iedb.org/>`_ to have your IP
   address removed from the IEDB blacklist.

**Support for local installation of the IEDB Analysis Resources**

pVACsplice provides the option of using a local installation of the IEDB MHC
`class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_
binding prediction tools.

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page (note: the pvactools 
   docker image now contains IEDB).

**MHC Class I and Class II predictions**

Both MHC Class I and Class II predictions are supported. Simply choose the desired
prediction algorithms and HLA alleles during processing and Class I and Class II
prediction results will be written to their own respective subdirectories in your
output directory. pVACsplice supports binding affinity algorithms as well as presentation
and immunogenicity algorithms.

By using the IEDB RESTful web interface, pVACsplice leverages their extensive support of different prediction algorithms.

In addition to IEDB-supported prediction algorithms, we've also added support
for a variety of additional algorithms.

.. list-table::
   :header-rows: 1

   * - Algorithm
     - Version(s)
     - MHC Class
     - Prediction Type
     - Supports Percentile Ranks?
     - Supports Normalized Percentile Ranks?
   * - BigMHC_EL
     -
     - MHC Class I
     - Presentation
     - no
     - yes
   * - BigMHC_IM
     -
     - MHC Class I
     - Immunogenicity
     - no
     - yes
   * - DeepImmuno
     -
     - MHC Class I
     - Immunogenicity
     - no
     - yes
   * - ImmuoScope_IM
     -
     - MHC Class II
     - Immunogenicity
     - no
     - no
   * - MHCflurry
     -
     - MHC Class I
     - Binding
     - yes
     - yes
   * - MHCflurryEL
     -
     - MHC Class I
     - Presentation, Processing
     - yes (Presentation only)
     - yes (Presentation and Processing)
   * - MHCnuggetsI
     -
     - MHC Class I
     - Binding
     - yes
     - yes
   * - MHCnuggetsII
     -
     - MHC Class II
     - Binding
     - yes
     - no
   * - MixMHC2pred
     -
     - MHC Class II
     - Presentation
     - yes
     - no
   * - MixMHCpred
     -
     - MHC Class I
     - Binding
     - yes
     - yes
   * - NNalign
     - 2.3
     - MHC Class II
     - Binding
     - yes
     - no
   * - NetMHC
     - 4.0
     - MHC Class I
     - Binding
     - yes
     - yes
   * - NetMHCIIpan
     - 4.0 (not supported by standalone IEDB), 4.1 (default), 4.2., 4.3
     - MHC Class II
     - Binding
     - yes
     - no
   * - NetMHCIIpanEL
     - 4.0 (not supported by standalone IEDB), 4.1 (default), 4.2., 4.3
     - MHC Class II
     - Presentation
     - yes
     - no
   * - NetMHCcons
     - 1.1
     - MHC Class I
     - Binding
     - yes
     - yes
   * - NetMHCpan
     - 4.1
     - MHC Class I
     - Binding
     - yes
     - yes
   * - NetMHCpanEL
     - 4.1
     - MHC Class I
     - Presentation
     - yes
     - yes
   * - PRIME
     -
     - MHC Class I
     - Immunogenicity
     - yes
     - yes
   * - Pickpocket
     - 1.1
     - MHC Class I
     - Binding
     - yes
     - yes
   * - SMM
     - 1.0
     - MHC Class I
     - Binding
     - yes
     - yes
   * - SMMPMBEC
     - 1.0
     - MHC Class I
     - Binding
     - yes
     - yes
   * - SMMalign
     - 1.1
     - MHC Class II
     - Binding
     - yes
     - no

**Calculation of normalized percentiles**

Not all prediction algorithms supported by pVACsplice output a percentile rank.
In order to alleviate this issue, and to provide percentile ranks that have been consistently
calculated, we have run predictions for all class I algorithms supported by pVACtools on 100,000
reference peptides each in lengths 8-11 and for the most common 1,000 human class I MHC alleles.
These predictions allow pVACsplice to support the calculation of normalized percentiles. This feature
is enable be setting the ``--use-normalized-percentiles`` parameter. With this option enabled,
pVACsplice will calculate normalized percentiles scores for all predicted neoantigen candidates and
selected prediction algorithms. These normalized percentile ranks will be used in place of percentile
ranks calculated by the algorithms natively. The algorithms supporting this feature are noted in the
table above.

**Comprehensive filtering**

Automatic filtering on the binding affinity IC50 (nm) value, binding percentile, presentation percentile,
and immunogenicity percentile narrows down the results to only include
"good" candidate peptides. The binding filter thresholds can be adjusted by the user for each
pVACsplice run. pVACsplice also support the option of filtering on allele-specific binding thresholds
as recommended by `IEDB <https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions>`_.
Additional filtering on the binding affitinity can be manually done by the user by running the
:ref:`standalone binding filter <pvacsplice_filter_commands>` on the filtered result file
to narrow down the candidate epitopes even further or on the unfiltered
all_epitopes file to apply different cutoffs.

Readcount and expression data are extracted from an annotated VCF to automatically filter with
adjustable thresholds on depth, VAF, and/or expression values. The user can also manually run
the :ref:`standalone coverage filter <pvacsplice_filter_commands>` to further narrow down their results
from the filtered output file.

If the input VCF is annotated with Ensembl transcript support levels (TSLs), MANE Select, and
Canonical status, pVACseq will filter on these to only keep high-confidence
transcripts. This filter can also be run :ref:`standalone
<pvacsplice_filter_commands>`.

As a last filtering step, pVACsplice applies the top score filter to only keep the top scoring epitope
for each variant. As with all previous filters, this filter can also be run
:ref:`standalone <pvacsplice_filter_commands>`. Please also see that section for more
details about how the top scoring epitope is determines.

**NetChop and NetMHCstab integration**

Cleavage position predictions are added with optional processing through NetChop.

Stability predictions can be added if desired by the user. These predictions are obtained via NetMHCstabpan.

**Reference proteome similarity analysis**

This optional feature will search for an epitope in the reference proteome
using BLAST or a reference proteome FASTA file to determine if the epitope occurs elsewhere in the proteome and
is, therefore, not tumor-specific.

**Problematic amino acids**

This optional feature allows users to specify a list of amino acids that would
be considered problematic to occur either everywhere or at specific positions
in a neoepitope. This can be useful when certain amino acids would be
problematic during peptide manufacturing.
