.. image:: ../../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _vep:

Annotating your VCF with VEP
============================

The input to the pVACseq pipeline is a VEP-annotated VCF. This will add
consequence, transcript, and gene information to your VCF.

Installing VEP
--------------

1. To download and install the VEP command line tool follow `the VEP installation instructions <http://useast.ensembl.org/info/docs/tools/vep/script/index.html>`_.
2. We recommend the use of the VEP cache for your annotation. The VEP cache
   can be downloaded following `these VEP cache installation instructions
   <http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache>`_.
   Please ensure that the Ensembl cache version matches the reference build
   and Ensembl version used in other parts of your analysis (e.g. for RNA-seq
   gene/transcript abundance estimation).
3. Download the VEP plugins from the `GitHub repository <https://github.com/Ensembl/VEP_plugins>`_
   by cloning the repository:

.. code-block:: none

   git clone https://github.com/Ensembl/VEP_plugins.git

4. :ref:`Copy the Wildtype and Frameshift plugins<install_vep_plugin_label>` provided with the
   pVACseq package to the folder with the other VEP plugins by running the following command:

.. code-block:: none

   pvacseq install_vep_plugin <VEP plugins directory>

Running VEP
-----------

**Example VEP Command**

.. code-block:: none

   ./vep \
   --input_file <input VCF> --output_file <output VCF> \
   --format vcf --vcf --symbol --terms SO --tsl --biotype \
   --hgvs --fasta <reference build FASTA file location> \
   --offline --cache [--dir_cache <VEP cache directory>] \
   --plugin Frameshift --plugin Wildtype \
   [--dir_plugins <VEP_plugins directory>] [--pick] [--transcript_version]

Required VEP Options
____________________

.. code-block:: none

   --format vcf
   --vcf
   --symbol
   --terms SO
   --tsl
   --biotype
   --hgvs
   --fasta <reference build FASTA location>
   --offline
   --cache
   --plugin Frameshift
   --plugin Wildtype

- The ``--format vcf`` option specifies that the input file is in VCF format.
- The ``--vcf`` option will result in the output being written in VCF format.
- The ``--symbol`` option will include gene symbol in the annotation.
- The ``--terms SO`` option will result in Sequence Ontology terms being used
  for the consequences.
- The ``--tsl`` option adds transcript support level information to the
  annotation.
- The ``--biotype`` option adds biotype of the transcript or regulatory
  feature to the annotation.
- The ``--hgvs`` option will result in HGVS identifiers being added to the
  annotation.
- Using the ``--hgvs`` option requires the usage of the ``--fasta`` argument to
  specify the location of the reference genome build FASTA file.
- The ``--offline`` option will eliminate all network connections for speed
  and/or privacy.
- The ``--cache`` option will result in the VEP cache being used for
  annotation.
- The ``--plugin Frameshift`` option will run the Frameshift plugin which will
  apply a frameshift mutation to a transcript sequence to compute the full mutated
  protein sequence.
- The ``--plugin Wildtype`` option will run the Wildtype plugin which will
  include the transcript protein sequence in the annotation.

Useful VEP Options
__________________

.. code-block:: none

   --dir_cache <VEP cache directory>
   --dir_plugins <VEP_plugins directory>
   --pick
   --transcript_version

- The ``--dir_cache <VEP cache directory>`` option may be needed if the VEP
  cache was downloaded to a different location than the default. The default
  location of the VEP cache is at ``$HOME/.vep``.
- The ``--dir_plugins <VEP_plugins directory>`` option may need to be set
  depending on where the VEP_plugins were installed to.
- The ``--pick`` option might be useful to limit the annotation to the "top"
  transcript for each variant (the one for which the most dramatic consequence 
  is predicted). Otherwise, VEP will annotate each variant with all possible
  transcripts. pVACseq will provide predictions for all transcripts in the VEP
  CSQ field. Running VEP without the ``--pick`` option can therefore drastically
  increase the runtime of pVACseq.
- The ``--transcript_version`` option will add the transcript version to the
  transcript identifiers. This option might be needed if you intend to
  annotate your VCF with expression information. Particularly if your
  expression estimation tool uses versioned transcript identifiers (e.g.
  ENST00000256474.2).

Additional VEP options that might be desired can be found
`here <http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html>`_.
