Installation
============

pVACviz is part of the pVACtools package. To install pVACtools, execute the following command on your Terminal:

.. code-block:: none

   pip install pvactools

More detailed installation instructions can be found :ref:`here <install>`. Note that the following are the bare minimum you 
need to run pVACviz. Most users will probably just want to complete the full pvactools install as described :ref:`here <install>`. 
That includes pVACviz along with all the other components, local installation of IEDB, etc.  You can also use the pvactools docker 
container which contains all tools and their dependencies (including those for pVACviz).

MHCflurry
---------

When installing pVACtools for the first time, you will need to manually
download the MHCflurry dataset:

.. code-block:: none

   mhcflurry-downloads fetch

PostgreSQL
----------

pVACviz requires a Postgres database. To install Postgres follow
the `installation instructions <http://postgresguide.com/setup/install.html>`_.

.. note::

   On Debian-based Linux distributions version Postgres V9.6 or lower is
   required.
