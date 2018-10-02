.. image:: ../images/pVACviz_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACviz logo

.. _pvacviz_visualizing:

Visualizing Processes
=====================

pVACviz provides a results visualization for exploring the results of pVACseq processes. It is able to visualize both the results from processes launched from pVACviz, and results from any pVACseq process.

Visualizing Completed Processes
-------------------------------
You may view visualizations of completed pVACseq processes launched from the pVACviz Start form from two locations within the application. 
The :ref:`Manage section <pvacviz_managing>` includes a process detail page, reachable by clicking on the Details link on the right side of rows in the process table.On the process detail page, the bottom right card contains a list of all files produced by the pVACseq process. Visualizable files will display a Visualize button that when clicked will load the visualization for that file.

Additionally, on the Visualize main page in the right hand column, all processes currently managed by pVACapi will be listed with their visualizable files. Clicking on a file will load the visualization for that file.

Visualizing pVACseq Results Files
---------------------------------
Any final results TSV file produced by pVACseq processes - whether launched via pVACviz or the command line - may be visualized with pVACviz. You may drop any file or folder in pVACapi's /visualize directory, and it will scan it for visualizable files. These files will then be listed on in the right column of the main Vizualize page. Click on any of the listed pages to launch the visualization

Scatterplot Visualization
-------------------------
On the visualization's scatterplot are placed all of the data points contained in the tsv results file, one dot per row. A set of tools along the right side of the visualization allow you to select and manipulate the plot in various ways.

These icons perform the following functions:

=================  ================  ========
Icon               Name              Function
=================  ================  ========
|pan|              **Pan**           The pan tool allows the user to pan the plot by left-dragging a mouse or dragging a finger across the plot region.
|box_zoom|         **Box Zoom**      The box zoom tool allows the user to define a rectangular region to zoom the plot bounds too, by left-dragging a mouse, or dragging a finger across the plot area.
|wheel_zoom|       **Wheel Zoom**    The wheel zoom tool will zoom the plot in and out, centered on the current mouse location. It will respect any min and max values and ranges preventing zooming in and out beyond these.
=================  ================  ========

.. |pan| image:: https://bokeh.pydata.org/en/latest/_images/Pan.png
   :align: middle
   :width: 16
   :height: 16

.. |box_zoom| image:: https://bokeh.pydata.org/en/latest/_images/BoxZoom.png
   :align: middle
   :width: 16
   :height: 16

.. |wheel_zoom| image:: https://bokeh.pydata.org/en/latest/_images/WheelZoom.png
   :align: middle
   :width: 16
   :height: 16

Filters
-------

Data Table
----------

Exporting Visualization Data
----------------------------

