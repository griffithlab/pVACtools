.. image:: ../images/pVACfuse_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACfuse logo

Prerequisites
=============

Fusion detection and annotation using AGFusion
----------------------------------------------

pVACfuse accepts a output directory from `AGFusion <https://github.com/murphycj/AGFusion>`_ as input.

AGFusion allows a user to annotate output files from several fusion callers
using the ``agfusion batch`` command. The below example is for annotating the
output from the STAR-Fusion caller but many other fusion callers are supported.
For a full list see the `AGFusion documentation <https://github.com/murphycj/AGFusion#input-from-fusion-finding-algorithms>`_.

.. code-block:: none

   agfusion batch \
   -f <star_fusion_tsv> \
   -a starfusion \
   -db agfusion.homo_sapiens.87.db \
   -o <output_directory> \
   --middlestar \
   --noncanonical

The ``--middlestar`` flag is required in order to use the ouput with pVACfuse.
This will indicate the fusion position in the fusion peptide sequence.

The ``--noncanonical`` flag is optional and can be used to annotate the fusion
with informations from all possible transcripts. By default only canonical
transcripts are used.

Fusion detection and annotation using Arriba
--------------------------------------------

Arriba detects fusions from STAR variant calls. For more information on how to
run STAR and Arriba see their `documentation <https://arriba.readthedocs.io/en/latest/workflow/#demo-script>`_.

Running STAR-Fusion for read support and expression information
---------------------------------------------------------------

In addition to predicting fusions, `STAR-Fusion <https://github.com/STAR-Fusion/STAR-Fusion>`_
will also provide read support
and expression information for the predicted fusions. The STAR-Fusion output
file (either ``star-fusion.fusion_predictions.tsv`` or
``star-fusion.fusion_predictions.abridged.tsv``) can be provided as an
optional input to pVACfuse using the ``--starfusion-file`` parameter.
