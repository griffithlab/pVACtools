.. image:: ../images/pVACview_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACview logo

.. raw:: html

  <style> .large {font-size: 70%; font-weight: bold} </style>
  <style> .bold {font-size: 110%; font-weight: bold} </style>
  <style> .underline {font-size: 100%; text-decoration: underline;} </style>

.. role:: large
.. role:: bold
.. role:: underline

.. _troubleshooting_pvacview_label:

Commonly Asked Questions
--------------------------

:bold:`What does each column of the main table mean?/ Where can I find detailed description of each column?`

Please reference to :ref:`pVACseq aggregated report output file documentation <aggregated>` for a description of each column.

These definitions as well as further details on the coloring and bar graphs in the main table can be found in the app by either: 1) hovering over individual column names and 2) clicking on the
tooltip at the top right of the main table.


:bold:`What does the peptide table show?`

Depending on the transcript you select, the peptide table will show the MT peptide sequences from that transcript that were predicted to be a good binder (along with their respective wildtype sequences). Currently, only predictions for good binding HLA alleles are displayed for MT/WT pairs and those with poor predictions are marked in X.

If the WT peptide sequence is unavailable in cases such as frameshift variants, these will be labelled with the mutation type followed by NA (e.g. FS-NA for frameshifts).


:bold:`What is the violin plot in the additional information box showing?`

pVACtools allows users to select up to 8 Class I algorithms and 4 Class II algorithms when making binding predictions. In the main table as well as the peptide table, the best peptide is determined using a median
value across all binding predictions from the selected algorithms. However `previous study <https://cancerimmunolres.aacrjournals.org/content/8/3/409>`_ has shown how binding algorithms can vary greatly in terms of its predictions for the same set of peptides. Thus, we want to
provide users with the ability to visualize the distribution of predictions across algorithms for each MT/WT pair.

The violin plots also provide a guiding line at 500nM and 1000nM for IC50 values and 0.5% and 2% for percentile values.

:bold:`What is the anchor plot in the additional information box showing?`

Anchor locations and its relative position can influence prioritization decisions for neoantigens as shown `here <https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1>`_ and we want to provide users
with the data to be able to take these considerations into account. In the anchor plot, peptide sequences (from the peptide table) are plotted and a heatmap overlays the sequences where a darker blue represents
a higher probability of being an anchor location. The mutation(s) is/are marked in red letters.

To the right of the additional information box, we provide a graphical guide regarding how one can take anchor information into account, more details on how to interpret this data can be found in `our paper <https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1>`_.

:bold:`I'm getting an error for the anchor heatmap tab saying "Error:polygon edge not found", what do I do?`

Users have occasionally ran into the problem where their anchor heatmap does not display and instead shows an error saying "polygon edge not found". After investigation
we believe this may be related to your arial font file. This stackoverflow page describes the detailed steps to resolving this issue:
`https://stackoverflow.com/questions/10581440/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon <https://stackoverflow.com/questions/10581440/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon>`_ .
