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

We recommend running pVACvector in various modes such as using median as well
as the lowest/best binding affinities across various epitope lengths against
the multiple prediction algorithms we offer. We also recommend setting a higher
binding threshold of 1000nM. By default, pVACvector runs on the median binding
affinity with a 500nM binding threshold. Running on the lowest/best binding affinity
ensures that if even just one of eight algorithms thinks a peptide is a good
binder, the junction will be rejected. However, running pVACvector with more
conservative parameters also increases the likelihood of pVACvector not being
able to find a valid path. Therefore, running pVACvector in these multiple modes,
and choosing the overall most conservative solution is recommended.

.. program-output:: pvacvector run -h

.. .. argparse::
        :module: tools.pvacvector.run
        :func: define_parser
        :prog: pvacvector run
