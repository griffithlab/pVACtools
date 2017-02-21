Features
========

**SNV and Indel support**

pVAC-Seq offers epitope binding predictions for missense, inframe indel, and frameshift mutations.

**VCF support**

pVAC-Seq uses a single-sample VCF file as its input. This VCF file must be annotated with VEP. See the :ref:`prerequisites_label` for more information.

**No local install of epitope prediction software needed**

pVAC-Seq utilizes the IEDB RESTful web interface. This means that none of the underlying prediction software, like NetMHC, needs to be installed locally.

.. warning::
   We only recommend using the RESTful API for small requests. If you use the
   RESTful API to process large VCFs or to make predictions for many alleles,
   epitope lengths, or prediction algorithms, you might overload their system.
   This can result in the blacklisting of your IP address by IEDB, causing
   403 errors when trying to use the RESTful API. In that case please open
   a ticket with `IEDB support <http://help.iedb.org/>`_ to have your IP
   address removed from the IEDB blacklist.

**Support for local installation of the IEDB Analysis Resources**

pVAC-Seq provides the option of using a local installation of the IEDB MHC
`class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_
binding prediction tools.

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

**MHC Class I and Class II predictions**

Both MHC Class I and Class II predictions are supported. Simply choose the desired prediction algorithms and HLA alleles during processing and Class I and Class II prediction results will be written to their own respective subdirectories in your output directory.

By using the IEDB RESTful web interface, pVAC-Seq leverages their extensive support of different prediction algorithms.

================================= =======
MHC Class I Prediction Algorithm  Version
================================= =======
NetMHCpan                         2.8
NetMHC                            4.0
NetMHCcons                        1.1
PickPocket                        1.1
SMM
SMMPMBEC
================================= =======

================================= =======
MHC Class II Prediction Algorithm Version
================================= =======
NetMHCIIpan                       3.0
SMMalign                          1.1
NNalign                           2.2
================================= =======

**Comprehensive filtering**

Automatic filtering on the binding affinity ic50 value narrows down the results to only include "good" candidate peptides. The binding filter threshold can be adjusted by the user for each pVAC-Seq run, and additional filtering can be manually done by the user on the final result file to narrow down the candidate epitopes even further.

bam-readcount and cufflinks files can be provided by the user as additional input files and are used to extract coverage and expression data. When any bam-readcount or cufflinks files are provided, automatic filtering with adjustable thresholds on depth, VAF, and/or expression value will narrow down the results. The user can also manually run the coverage filter to further narrow down their results from the final output file.

The user can also specify an option to only keep the top scoring result for each allele-peptide length combination for each variant.

**NetChop and NetMHCstab integration**

Cleavage position predictions are added with optional processing through NetChop.

Stability predictions can be added if desired by the user. These predictions are obtained via NetMHCstab.
