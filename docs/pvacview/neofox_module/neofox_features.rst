.. image:: ../../images/pVACview_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACview logo

.. _neofox_features:

.. raw:: html

  <style> .large {font-size: 90%; font-weight: bold} </style>
  <style> .bold {font-size: 100%; font-weight: bold} </style>

.. role:: large
.. role:: bold

Neofox Features
---------------

:large:`Additional Modules for Viewing Data`
_________________________________________________________________________________________________________________________

To expand the usage of pVACview we include two additional modules for viewing neoantigen prediction data: Neofox and Custom.

- :bold:`Neofox:`

  NeoFox is a Python package that annotates neoantigen candidates with published neoantigen features (citation??). 
  It is a popular neoantigen evaluation tool and this module enables a user to flexibly visualize NeoFox output data.
  The main table displays all annotated neoantigens and selecting a candidate results in violin plots to be 
  generated which shows the selected candidate in red compared to all other candidates. Up to 6 features can be viewed as violin plots at a time.

  .. figure:: ../../images/screenshots/pvacview-neofox-table-violinplot.png
              :width: 1000px
              :align: left
              :figclass: align-left

  You can also further investigate the data using the dynamic scatter plot where you can choose any feature to be the X-axis, Y-axis,
  Color, or Size variable. The X and Y scale can be transformed and a range of values subsetted. The color represents the minimum
  and maximum values can also be changed to any HEX value. 

  .. figure:: ../../images/screenshots/pvacview-neofox-dynamic-scatter.png
              :width: 1000px
              :align: left
              :figclass: align-left

  We have also marked features our lab finds particularly interesting with an asterisk. These are features that we would place high value on
  when selecting the best neoantigens. 

