.. image:: images/pVACviz_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACviz logo

.. _pvacviz:

pVACviz
=======

pVACviz is the browser-based, graphical user interface for the pVACtools command line tools.
It currently supports starting and managing :ref:`pvacseq` runs as well
as visualizing the results of your runs.

Running pVACviz
---------------

To run pVACviz you first need to start pVACapi, which is used to communicate
between the user interface and the command line tool. pVACapi can be started
by executing the following command on the command line:

.. code::

   pvacapi

Depending on the number of completed and running processes it must check, the API may take several seconds to start up. After pVACapi has started, launch pVACviz in a separate terminal window by executing the following command on the command line:

.. code::

   pvacviz

The command will start a HTTP server that provides the web client files and assets, and opens up the client in the default web browser specified by your operating system. In some cases, pvacviz will not be able to automatically open the web browser. If no browser launches after starting pvacviz, you will need to manually load the URL, `http://localhost:4200/`, in a Firefox, Chrome, or Safari browser.

.. toctree::
   :maxdepth: 2
   :glob:

   pvacviz/directories
   pvacviz/starting
   pvacviz/managing
   pvacviz/visualizing
