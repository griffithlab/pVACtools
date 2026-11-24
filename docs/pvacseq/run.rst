.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. _pvacseq_run:

Usage
====================================

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

.. program-output:: pvacseq run -h

ML Accept Preselection
----------------------

With ``--run-ml-predictions``, at most 20 ML-predicted Accept candidates per
sample are preselected as ``Accept`` in the ``Evaluation`` column by default.
Candidates are ranked by their full-precision Accept probability. Additional
predicted Accept candidates are set to ``Pending``.

Use ``--ml-accept-limit N`` to set a numeric limit or
``--ml-accept-limit none`` to preselect every predicted Accept candidate.
With ``--ml-accept-limit 0``, all predicted Accept candidates are set to
``Pending`` in ``Evaluation``. This option also applies to
``pvacseq add_ml_predictions``.

..  .. argparse::
        :module: tools.pvacseq.run
        :func: define_parser
        :prog: pvacseq run
