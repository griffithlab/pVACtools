.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

Annotating your VCF with VEP
============================

The input to the pVACseq pipeline is a VEP-annotated VCF. This will add
consequence, transcript, and gene information to your VCF.

Installing VEP
--------------

1. To download and install the VEP command line tool follow `these instructions <http://useast.ensembl.org/info/docs/tools/vep/script/index.html>`_.
2. Download the VEP_plugins from their `GitHub repository <https://github.com/Ensembl/VEP_plugins>`_.
3. :ref:`Copy the Wildtype plugin<install_vep_plugin_label>` provided with the pVACseq package to the folder with the other VEP_plugins:

.. code-block:: none

   pvacseq install_vep_plugin

Running VEP
-----------

**Example VEP Command**

.. code-block:: none

   perl variant_effect_predictor.pl \
   --input_file <input VCF> --format vcf --output_file <output VCF> \
   --vcf --symbol --terms SO --hgvs --plugin Downstream --plugin Wildtype \
   [--dir_plugins <VEP_plugins directory>] [--pick] [--transcript_version]

Required VEP Options
____________________

.. code-block:: none

   --format vcf
   --vcf
   --symbol
   --plugin Downstream
   --plugin Wildtype
   --terms SO
   --hgsv

- The ``--format vcf`` option specifies that the input file is in VCF format.
- The ``--vcf`` option will result in the output being written in VCF format.
- The ``--symbol`` option will include gene symbol in the annotation.
- The ``--plugin Downstream`` option will run the Downstream plugin which will
  compute the downstream protein sequence after a frameshift.
- The ``--plugin Wildtype`` option will run the Wildtype plugin which will
  include the transcript protein sequence in the annotation.
- The ``--terms SO`` option will result in Sequence Ontology terms being used
  for the consequences.
- The ``--hgvs`` option will result in HGVS identifiers being added to the
  annotation.

Useful VEP Options
__________________

.. code-block:: none

   --dir_plugins <VEP_plugins directory>
   --pick
   --transcript_version

- The ``--dir_plugins <VEP_plugins directory>`` option may need to be set
  depending on where the VEP_plugins were installed to.
- The ``--pick`` option might be useful to limit the annotation to the top
  transcripts. Otherwise, VEP will annotate each variant with all possible
  transcripts. pVACseq will provide predictions for all transcripts in the VEP
  CSQ field. Running VEP without the ``--pick`` option can therefor drastically
  increase the runtime of pVACseq.
- The ``--transcript_version`` option will add the transcript version to the
  transcript identifiers. This option might be needed when annotation the VCF
  further with expression information and the expression tool uses versioned
  transcripts in their identifier.

Additional VEP options that might be desired can be found
`here <http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html>`_.
