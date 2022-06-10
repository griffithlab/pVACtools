.. image:: ../images/pVACview_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACview logo

.. _pvacview_prerequisites:

.. raw:: html

  <style> .large {font-size: 90%; font-weight: bold} </style>
  <style> .bold {font-size: 100%; font-weight: bold} </style>

.. role:: large
.. role:: bold

Prerequisites
---------------

In order to launch the pVACview R shiny application, you will need to have R/ R studio and a list of R packages correctly installed.
Once launched, pVACview will require you to upload your corresponding input files for analysis. Alternatively, there is also a demo dataset that can be loaded to explore basic features of the app.
An online version of pVACview is also available at `pvacview.org <https://www.pvacview.org>`_.

:bold:`Please note that you will need internet connection for pVACview, as it needs to download necessary datasets including both the demo data and anchor calculations.`

Installing R/ R studio
____________________________

In order to use pVACview, you will need to download R. Please refer `here <https://cran.rstudio.com/>`_ for downloading R (version 3.5 and above required).
You may also take the additional step of `downloading R studio <https://www.rstudio.com/products/rstudio/download/>`_ if you are not familiar with launching R Shiny from the command line.

Additionally, there are a number of packages you will need to install in your R/R studio, instructions are provided below:

.. code-block:: none

  install.packages("shiny", dependencies=TRUE)
  install.packages("ggplot2", dependencies=TRUE)
  install.packages("DT", dependencies=TRUE)
  install.packages("reshape2", dependencies=TRUE)
  install.packages("jsonlite", dependencies=TRUE)
  install.packages("tibble", dependencies=TRUE)
  install.packages("tidyr", dependencies=TRUE)
  install.packages("plyr", dependencies=TRUE)
  install.packages("dplyr", dependencies=TRUE)
  install.packages("shinydashboard", dependencies=TRUE)
  install.packages("shinydashboardPlus", dependencies=TRUE)
  install.packages("fresh", dependencies=TRUE)
  install.packages("shinycssloaders", dependencies=TRUE)
  install.packages("RCurl", dependencies=TRUE)
  install.packages("curl", dependencies=TRUE)


Note that certain R packages may have further dependencies that require additional installation. If you are using R studio, these should be automatically installed for you during the respective package
installation process. However, if you are using conda-based R, you may need to install them manually (usually by running ``install.packages(<package name>)``).

.. _launching_pvacview_label:

Launching pVACview R Shiny App
______________________________

Once you have R packages and their respective dependencies successfully installed, you may now launch the pVACview app.

.. role:: large

:large:`Option 1:`

If using R from the command line, you will first need to locate the path to the shiny app e.g. ``~/pVACtools_output/MHC_Class_I/``. Then you can run the following command from your
terminal or console window once you replace ``~/pVACtools_output/MHC_Class_I/`` with the appropriate path:

.. code-block:: none

  R -e "shiny::runApp('~/pVACtools_output/MHC_Class_I/')"

Your terminal/console window should then provide you with a localhost port where the R shiny app is being hosted. e.g. ``http://127.0.0.1:6126`` or ``http://localhost:6126``. Open a browser window of
your choice and enter the address to use pVACview.

:large:`Option 2:`

Alternatively, you can run R studio and open the file ``app.R`` located in your pVACtools output folder. In your R studio interface, on the top right, you should see a ``runApp`` button and
upon clicking pVACview will be launched. You may also run the following line in R studio console to achieve the same result:

.. code-block:: none

  runApp()

This will automatically launch a browser window from R studio, where the local port is specified on the top left corner (you can also find it in the R studio console e.g. ``Listening on http://127.0.0.1:6126``).
To ensure full functionality, you should take the local port address and launch it in a browser of your choice (e.g. Chrome, firefox etc).
