Installation
============

pVACtools is written for Linux but some users have been able to run it successfully on Mac OS X. If you are using Windows you will need to set up a Linux environment, for example by setting up a virtual machine.

pVACtools requires Python 3.5. Before running any installation steps, check the Python version installed on your system:

.. code-block:: none

   python -V

If you don't have Python 3.5 installed, we recommend using `Conda <http://conda.pydata.org/docs/py2or3.html>`_ to emulate a Python 3.5. environment. We've encountered problems with users that already have Python 2.x installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.x installed we **strongly** recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use ``pip`` to install pVACtools. Make sure you have ``pip`` installed. ``pip`` is generally included in python distributions, but may need to be upgraded before use. See the `instructions <https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel>`_ for installing or upgrading ``pip``.

After you have pip installed, type the following command on your Terminal:

.. code-block:: none

   pip install pvactools

You can check that ``pvactools`` has been installed under the default environment by listing all installed packages:

.. code-block:: none

   pip list

``pip`` will fetch and install pVACtools and its dependencies for you. After installing, each tool of the pVACtools suite is available in its own command line tree directly from the Terminal.

If you have an old version of pVACtools installed you might want to consider upgrading to the latest version:

.. code-block:: none

   pip install pvactools --upgrade

.. _iedb_install:

Installing IEDB binding prediction tools (strongly recommended)
---------------------------------------------------------------

.. warning::
   Using a local IEDB installation is strongly recommended for larger datasets
   or when the making predictions for many alleles, epitope lengths, or
   prediction algorithms.

.. warning::
   The IEDB binding prediction tools are only compatible with Linux.

You may create a local install of the IEDB binding prediction tools by first downloading the archives for `class I <http://tools.iedb.org/mhci/download/>`_ and `class II <http://tools.iedb.org/mhcii/download/>`_ from the IEDB website. If using both the Class I and the Class II tools, they both need to be installed into the same parent directory.
   
.. note::

   IEDB requires tcsh. You can install it by running ``sudo apt-get install tcsh``.

MHC Class I
___________

Download the archives for `class I <http://tools.iedb.org/mhci/download/>`_ and unpack them.

.. code-block:: none

   tar -zxvf IEDB_MHC_I-2.17.tar.gz
   cd mhc_i
   ./configure
    
.. note::

   Running the ``configure`` script requires a Python 2 environment. If you are currently emulating a Python 3 environment with Conda you will need to run ``source deactivate`` before executing the ``configure`` script.

Open ``method/netmhc-4.0-executable/netmhc_4_0_executable/__init__.py`` and delete/comment out the first line (``import pkg_resources``). Also delete/comment out the same line of code from ``method/netmhcpan-3.0-executable/netmhcpan_3_0_executable/__init__.py`` on line 7.

If you want to use the NetMHCcons prediction algorithm you will need to change the shebang line of certain files to explicitly use python2.7. The files in question are:

* ``method/netmhccons-1.1-executable/netmhccons_1_1_executable/bin/pseudofind``
* ``method/netmhc-3.4-executable/netmhc_3_4_executable/netMHC``

In these files change the shebang line to ``#! /usr/bin/env python2.7``.

MHC Class II
____________

Download the archives for `class II <http://tools.iedb.org/mhcii/download/>`_ and unpack them.

.. code-block:: none

   tar -zxvf IEDB_MHC_II-2.17.3.tar.gz
   cd mhc_ii

On older versions of the IEDB software, you might need to update some paths in the configure scripts to use relative paths. Open the ``configure.py`` file and update the lines that set the ``smm`` and ``nn`` variables to use relative paths like so:

.. code-block:: none

   smm = re.compile(curDir + "/netMHCII-1.1")
   nn = re.compile(curDir + "/netMHCII-2.2")

Then run the configure script.

.. code-block:: none

   ./configure.py

.. note::

   Running the ``configure`` script requires a Python 2 environment. If you are currently emulating a Python 3 environment with Conda you will need to run ``source deactivate`` before executing the ``configure`` script.

Docker and CWL
--------------

A Docker container for pVACtools is available on DockerHub using the
`mgibio/pvactools <https://hub.docker.com/r/mgibio/pvactools/>`_ repo. CWL
tool wrappers for pVACseq, pVACfuse, and pVACvector can be downloaded
using the ``pvactools download_cwls`` command.
