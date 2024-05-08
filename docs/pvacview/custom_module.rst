.. image:: ../images/pVACview_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACview logo

.. _custom_module:

Custom Module
---------------

The custom module boasts the most flexibility for viewing neoantigen candidate data, since there are no required features that
are expected to be in the file.

We provide three examples of neoantigen prediction pipeline output data:

**Vaxrank**

Therapeutic vaccines targeting mutant tumor antigens (“neoantigens”) are an increasingly popular form of 
personalized cancer immunotherapy. Vaxrank is a computational tool for selecting neoantigen vaccine peptides 
from tumor mutations, tumor RNA data, and patient HLA type. Vaxrank is freely available on `github
<www.github.com/openvax/vaxrank>`_ under the Apache 2.0 open source license and 
can also be installed from the Python Package Index.

**NeoPredPipe**

`NeoPredPipe (Neoantigen Prediction Pipeline) <https://github.com/MathOnco/NeoPredPipe>`_
is offered as a contiguous means of predicting putative 
neoantigens and their corresponding recognition potentials for both single and multi-region tumor samples. 
This tool allows a user to process neoantigens predicted from single- or multi-region vcf files using ANNOVAR 
and netMHCpan.

**antigen.garnish.2**

Human and mouse ensemble tumor neoantigen prediction from SNVs and complex variants. 
Immunogenicity filtering based on the Tumor Neoantigen Selection Alliance (TESLA).
https://github.com/andrewrech/antigen.garnish

.. toctree::
   :maxdepth: 2
   :glob:

   custom_module/custom_upload
   custom_module/custom_features
   
