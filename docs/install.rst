Installation
============

pVAC-Seq is written for Linux but some users have been able to run it successfully on Mac OS X. If you are using Windows you will need to set up a Linux environment, for example by setting up a virtual machine.

pVAC-Seq requires Python 3.5. Before running any installation steps check the Python version installed on your system:

.. code-block:: none

   python -V

If you don't have Python 3.5 installed, we recommend using `Conda <http://conda.pydata.org/docs/py2or3.html>`_ to emulate a Python 3.5. environment. We've encountered problems with users that already have Python 2.x installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.x installed we **strongly** recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use ``pip`` to install pVAC-Seq. Make sure you have ``pip`` installed. ``pip`` is generally included in python distributions, but may need to be upgraded before use. See the `instructions <https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel>`_ for installing or upgrading ``pip``.

After you have pip installed, type the following command on your Terminal (for Mac and Linux users) or the Command Prompt (for Windows users):

.. code-block:: none

   pip install pvacseq

You can check that ``pvacseq`` has been installed under the default environment by listing all installed packages:

.. code-block:: none

   pip list

``pip`` will fetch and install pVAC-Seq and its dependencies for you. After installing, you can run ``pvacseq`` directly from the Terminal/Command Prompt.

If you have an old version of pVAC-Seq installed you might want to consider upgrading to the latest version:

.. code-block:: none

   pip install pvacseq --upgrade

.. _iedb_install:

Installing IEDB binding prediction tools (strongly recommended)
---------------------------------------------------------------

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms. More information on how to install IEDB locally can
   be found on the :ref:`Installation <iedb_install>` page.

You may create a local install of the IEDB binding prediction tools by first downloading the archives for `class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_ from the IEDB website. If using both the Class I and the Class II tools, they both need to be installed into the same parent directory.
   
.. note::

   IEDB requires tcsh. You can install it by running ``sudo apt-get install tcsh``.

MHC Class I
___________

Download the archives for `class I <http://tools.iedb.org/mhci/download/>`_ and unpack them.

.. code-block:: none

   tar -zxvf IEDB_MHC_I-2.15.2.tar.gz
   cd mhc_i
   ./configure
    
.. note::

   Running the ``configure`` script requires a Python 2 environment. If you are currently emulating a Python 3 environment with Conda you will need to run ``source deactivate`` before executing the ``configure`` script.

Open ``method/netmhc_4_0_executable/__init__.py`` and delete/comment out the first line (``import pkg_resources``). Also delete/comment out the same line of code from ``method/netmhcpan_3_0_executable/__init__.py`` on line 7.

If you want to use the NetMHCcons prediction algorithm you will need to change the shebang line of certain files to explicitly use python2.7. The files in question are:

* ``method/netMHCcons-1.1/bin/pseudofind``
* ``method/netMHC-3.4/netMHC``

In these files change the shebang line to ``#! /usr/bin/env python2.7``.

MHC Class II
____________

.. code-block:: none

   tar -zxvf IEDB_MHC_II-2.16.tar.gz
   cd mhc_ii
   ./configure.py
    
.. note::

   Running the ``configure`` script requires a Python 2 environment. If you are currently emulating a Python 3 environment with Conda you will need to run ``source deactivate`` before executing the ``configure`` script.
