.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Features
========

**SNV and Indel support**


pVACseq offers epitope binding predictions for missense, in-frame insertion, in-frame deletion, protein-altering, and frameshift mutations.

**VCF support**

pVACseq uses a VCF file as its input. This VCF file must contain sample genotype information and be annotated with the Ensembl Variant Effect Predictor (VEP). See the :ref:`prerequisites_label` section for more information.

**No local install of epitope prediction software needed**

pVACseq utilizes the IEDB RESTful web interface. This means that none of the underlying prediction software, like NetMHC, needs to be installed locally.

.. warning::
   We only recommend using the RESTful API for small requests. If you use the
   RESTful API to process large VCFs or to make predictions for many alleles,
   epitope lengths, or prediction algorithms, you might overload their system.
   This can result in the blacklisting of your IP address by IEDB, causing
   403 errors when trying to use the RESTful API. In that case please open
   a ticket with `IEDB support <http://help.iedb.org/>`_ to have your IP
   address removed from the IEDB blacklist.

**Support for local installation of the IEDB Analysis Resources**

pVACseq provides the option of using a local installation of the IEDB MHC
`class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_
binding prediction tools.

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page (note: the pvactools 
   docker image now contains IEDB).

**MHC Class I and Class II predictions**

Both MHC Class I and Class II predictions are supported. Simply choose the desired prediction algorithms and HLA alleles during processing and Class I and Class II prediction results will be written to their own respective subdirectories in your output directory.

By using the IEDB RESTful web interface, pVACseq leverages their extensive support of different prediction algorithms.

In addition to IEDB-supported prediction algorithms, we've also added support
for `MHCflurry <http://www.biorxiv.org/content/early/2017/08/09/174243>`_ and
`MHCnuggets <http://karchinlab.org/apps/appMHCnuggets.html>`_.

================================= =======
MHC Class I Prediction Algorithm  Version
================================= =======
NetMHCpan                         2.8
NetMHC                            4.0
NetMHCcons                        1.1
PickPocket                        1.1
SMM
SMMPMBEC
MHCflurry
MHCnuggets
================================= =======

================================= =======
MHC Class II Prediction Algorithm Version
================================= =======
NetMHCIIpan                       3.0
SMMalign                          1.1
NNalign                           2.2
MHCnuggets
================================= =======

**Comprehensive filtering**

Automatic filtering on the binding affinity ic50 (nm) value narrows down the results to only include
"good" candidate peptides. The binding filter threshold can be adjusted by the user for each
pVACseq run. pVACseq also support the option of filtering on allele-specific binding thresholds
as recommended by `IEDB <https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions>`_.
Additional filtering on the binding affitinity can be manually done by the user by running the
:ref:`standalone binding filter <filter_commands>` on the filtered result file
to narrow down the candidate epitopes even further or on the unfiltered
all_epitopes file to apply different cutoffs.

Readcount and expression data are extracted from an annotated VCF to automatically filter with
adjustable thresholds on depth, VAF, and/or expression values. The user can also manually run
the :ref:`standalone coverage filter <filter_commands>` to further narrow down their results
from the filtered output file.

If the input VCF is annotated with Ensembl transcript support levels (TSLs), pVACseq will
filter on the transcript support level to only keep high-confidence
transcripts of level 1. This filter can also be run :ref:`standalone
<filter_commands>`.

As a last filtering step, pVACseq applies the top score filter to only keep the top scoring epitope
for each variant. As with all previous filters, this filter can also be run
:ref:`standalone <filter_commands>`.

**Scoring of candidate neoepitopes**

Filtered neoepitopes are :ref:`scored and ranked <score>` based on the binding affinity,
fold change between mutant and wildtype binding affinity (agretopicity), gene expression, RNA
and DNA VAF.

**Incorporation of proximal germline and somatic variants**

To incorporate proximal variants into the neoepitope predictions, users can provide 
a :ref:`phased VCF of proximal variants <phased_vcf>` as an input to their pVACseq runs. 
This VCF is then used to incorporate amino acid changes of nearby variants that are in-phase 
with a somatic variant of interest. This results in corrected mutant and wildtype 
protein sequences that account for proximal variants when MHC binding predictions are performed.

**NetChop and NetMHCstab integration**

Cleavage position predictions are added with optional processing through NetChop.

Stability predictions can be added if desired by the user. These predictions are obtained via NetMHCstabpan.
