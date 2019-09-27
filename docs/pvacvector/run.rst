.. image:: ../images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

.. _pvacvector_run:

Usage
====================================

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

It may be necessary to explore the parameter space a bit when running pVACvector.
As binding predictions for some sites vary substantially across algorithms, the
most conservative settings may result in no valid paths, often due to one
"outlier" prediction. Carefully choosing which predictors to run may help
ameliorate this issue as well.

In general, setting a lower binding threshold (e.g., 500nM) and using the median
binding value (``--top-score-metric median``) will lead to greater possibility
of a design, while more conservative settings of 1000nM and lowest/best binding
value (``--top-score-metric lowest``) will give more confidence that there are
no junctional neoepitopes.

Our current recommendation is to run pVACvector several different ways, and
choose the path resulting from the most conservative set of parameters.

.. program-output:: pvacvector run -h

.. .. argparse::
        :module: tools.pvacvector.run
        :func: define_parser
        :prog: pvacvector run
