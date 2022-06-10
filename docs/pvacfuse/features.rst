.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Features
========

**Fusion analysis**

pVACfuses proceses fusion variants annotated by AGFusion.

**No local install of epitope prediction software needed**

pVACbind utilizes the IEDB RESTful web interface. This means that none of the underlying prediction software, like NetMHC, needs to be installed locally.

.. warning::
   We only recommend using the RESTful API for small requests. If you use the
   RESTful API to process large VCFs or to make predictions for many alleles,
   epitope lengths, or prediction algorithms, you might overload their system.
   This can result in the blacklisting of your IP address by IEDB, causing
   403 errors when trying to use the RESTful API. In that case please open
   a ticket with `IEDB support <http://help.iedb.org/>`_ to have your IP
   address removed from the IEDB blacklist.

**Support for local installation of the IEDB Analysis Resources**

pVACbind provides the option of using a local installation of the IEDB MHC
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
NetMHCpan                         4.1
NetMHC                            4.0
NetMHCcons                        1.1
PickPocket                        1.1
SMM                               1.0
SMMPMBEC                          1.0
MHCflurry
MHCnuggets
================================= =======

================================= =======
MHC Class II Prediction Algorithm Version
================================= =======
NetMHCIIpan                       4.0
SMMalign                          1.1
NNalign                           2.3
MHCnuggets
================================= =======

**Comprehensive filtering**

Automatic filtering on the binding affinity ic50 (nm) value narrows down the results to only include
"good" candidate peptides. The binding filter threshold can be adjusted by the user for each
pVACfuse run. pVACfuse also support the option of filtering on allele-specific binding thresholds
as recommended by `IEDB <https://help.iedb.org/hc/en-us/articles/114094151811-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions>`_.
Additional filtering on the binding affitinity can be manually done by the user by running the
:ref:`standalone binding filter <filter_commands>` on the filtered result file
to narrow down the candidate epitopes even further or on the unfiltered
all_epitopes file to apply different cutoffs.

pVACfuse also runs a top score filter to only keep the top scoring epitope
for each fusion. This filter can also be run
:ref:`standalone <filter_commands>`.

**NetChop and NetMHCstab integration**

Cleavage position predictions are added with optional processing through NetChop.

Stability predictions can be added if desired by the user. These predictions are obtained via NetMHCstabpan.

**Reference proteome similarity analysis**

This optional feature will search for an epitope in the reference proteome
using BLAST to determine if the epitope occurs elsewhere in the proteome and
is, therefore, not tumor-specific.
