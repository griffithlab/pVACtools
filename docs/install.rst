.. _install:

Installation
============

pVACtools is written for Linux but some users have been able to run it successfully on Mac OS X. If you are using Windows you will need to set up a Linux environment, for example by setting up a virtual machine.

pVACtools requires Python 3.6 or above. Before running any installation steps, check the Python version installed on your system:

.. code-block:: none

   python -V

If you don't have Python 3 installed, we recommend using `Conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ to emulate a Python 3 environment. We've encountered problems with users that already have Python 2.x installed when they also try to install Python 3. The defaults will not be set correctly in that case. If you already have Python 2.x installed we **strongly** recommmend using Conda instead of installing Python 3 locally.

Once you have set up your Python 3 environment correctly you can use ``pip`` to install pVACtools. Make sure you have ``pip`` installed. ``pip`` is generally included in python distributions, but may need to be upgraded before use. See the `instructions <https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel>`_ for installing or upgrading ``pip``.

After you have pip installed, type the following command on your Terminal:

.. code-block:: none

   pip install pvactools

You can check that ``pvactools`` has been installed under the default environment like so:

.. code-block:: none

   pip show pvactools

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

pVACtools is only compatible with IEDB 3.1 and above. We have tested pVACtools with the versions of IEDB class I and II listed below. Using a different version may cause errors.

.. important::
   By using the IEDB software, you are consenting to be bound by and become a
   "Licensee" for the use of IEDB tools and are consenting to the terms and
   conditions of the Non-Profit Open Software License ("Non-Profit OSL") version 3.0.

   Please read these two license agreements `here <http://tools.iedb.org/mhci/download/>`_
   before proceeding. If you do not agree to all of the terms of these two agreements,
   you must not install or use the product. Companies (for-profit entities) interested
   in downloading the command-line versions of the IEDB tools or running the entire analysis
   resource locally, should contact IEDB (license@iedb.org) for details on licensing options.

   Citing the IEDB

   All publications or presentations of data generated by use of the IEDB
   Resource Analysis tools should include citations to the relevant reference(s),
   found `here <http://tools.iedb.org/mhci/reference/>`_.

MHC Class I
___________

Download the archives for `class I <http://tools.iedb.org/mhci/download/>`_ and unpack them.

.. code-block:: none

   apt-get update && apt-get install -y tcsh gawk
   wget https://downloads.iedb.org/tools/mhci/3.1.5/IEDB_MHC_I-3.1.5.tar.gz
   tar -zxvf IEDB_MHC_I-3.1.5.tar.gz
   cd mhc_i
   ./configure

MHC Class II
____________

Download the archives for `class II <http://tools.iedb.org/mhcii/download/>`_ and unpack them.

.. code-block:: none

   apt-get update && apt-get install -y tcsh gawk
   wget https://downloads.iedb.org/tools/mhcii/3.1.12/IEDB_MHC_II-3.1.12.tar.gz
   tar -zxvf IEDB_MHC_II-3.1.12.tar.gz
   cd mhc_ii
   ./configure.py


Then run the configure script.

.. code-block:: none

   ./configure.py


Installing MHCflurry
--------------------

If you wish to run the MHCflurry prediction algorithm, you will need to
install the ``mhcflurry`` python package on your system. This package is set
as a dependency for the ``pvactools`` package so it should be installed
automatically when you download or upgrade the ``pvactools`` package. You can
install it manually by running:

.. code-block:: none

   pip install mhcflurry

.. note::

   The ``mhcflurry`` package needs to be installed in the same Python 3 environment as the ``pvactools`` package.

Next, you will need to download the download the MHCflurry datasets and trained models:

.. code-block:: none

   mhcflurry-downloads fetch

.. note::

   The ``mhcflurry-downloads fetch`` command will need to be run manually, even
   if the mhcflurry package was already installed automatically as a
   dependency with the ``pvactools`` package.

You can check that the ``mhcflurry`` package was installed successfully by running:

.. code-block:: none

  mhcflurry-predict -h

This should pull up the help page for the MHCflurry predictor.

Please note that MHCflurry depends on tensorflow, which will automatically be installed as a
dependency to the ``mhcflurry`` package. Newer versions of tensorflow might not be compatible
with older CPUs. In that case you will see a core dump failure. Downgrading
tensorflow manually to version 1.5.0 should solve this problem:

.. code-block:: none

   pip install tensorflow==1.5.0

Installing MHCnuggets
---------------------

If you wish to run the MHCnuggets prediction algorithm, you will need to
install the ``mhcnuggets`` python package on your system. This package is set
as a dependency for the ``pvactools`` package so it should be installed
automatically when you download or upgrade the ``pvactools`` package. You can
install it manually by running:

.. code-block:: none

   pip install mhcnuggets

.. note::

   The ``mhcnuggets`` package needs to be installed in the same Python 3
   environment as the ``pvactools`` package.

You can check that the ``mhcnuggets`` package was installed successfully by running:

.. code-block:: none

   pip show mhcnuggets

This should show information about the mhcnuggets package.

Please note that MHCnuggets depends on tensorflow, which will automatically be installed as a
dependency to the ``mhcnuggets`` package. Newer versions of tensorflow might not be compatible
with older CPUs. In that case you will see a core dump failure. Downgrading
tensorflow manually to version 1.5.0 should solve this problem:

.. code-block:: none

   pip install tensorflow==1.5.0

Installing BigMHC
-----------------

If you wish to run the BigMHC_EL or BigMHC_IM prediction algorithms, you will need to
install BigMHC on your system. This package not a direct dependency of
the the ``pvactools`` packages and needs to be installed manually by running:

.. code-block:: none

   pip install git+https://github.com/griffithlab/bigmhc.git#egg=bigmhc

.. note::

   BigMHC needs to be installed in the same python 3
   environment as the ``pvactools`` package.

You can check that BigMHC was installed successfully by running:

.. code-block:: none

   pip show bigmhc

This should show information about the BigMHC installation.

Installing DeepImmuno
---------------------

If you wish to run the DeepImmuno prediction algorithm, you will need to
install DeepImmuno on your system. This package not a direct dependency of
the the ``pvactools`` packages and needs to be installed manually by running:

.. code-block:: none

   pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno

.. note::

   DeepImmuno needs to be installed in the same python 3
   environment as the ``pvactools`` package.

You can check that DeepImmuno was installed successfully by running:

.. code-block:: none

   pip show deepimmuno

This should show information about the DeepImmuno installation.

.. _blast:

Installing BLAST
----------------

To run the reference proteome similarity step, standalone BLAST may be used.
To install BLAST please see `the official documentation
<https://www.ncbi.nlm.nih.gov/books/NBK52640/>`_. The BLAST tool needed is Protein BLAST
(blastp). Please make note of the installation path of blastp (retrievable by
calling ``which blastp``), as that path is needed for the ``--blastp-path`` argument in
the various pVACtools commands.

You will also need to install either the ``refseq_select_prot`` or the
``refseq_protein`` BLAST reference proteome databases. You can do so by running the
``update_blastdb.pl`` script provided with your BLAST installation (located in
the ``bin`` subdirectory). You will need to set the ``BLASTDB`` to point to the
installation directory of your BLAST reference proteome databases.

Downloading Reference Proteome FASTA file
-----------------------------------------

As an alternative to BLAST, a reference proteome fasta file may be used for
the reference proteome similarity step and specified as an input via the
``--peptide-fasta`` command. Any proteome fasta may be used. Ensembl provides
reference proteome fastas for many species. For example, the latest reference
proteome fasta for human can be downloaded from `this
link <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz>`_.

Docker and CWL/WDL
------------------

Versioned Docker containers for pVACtools are available on DockerHub using the
`griffithlab/pvactools <https://hub.docker.com/r/griffithlab/pvactools/>`_ repo.
The Docker container contains pVACtools as well as installations of the
standalone IEDB MHC Class I and Class II software and all other supported prediction
algorithms (MHCflurry, MHCnuggets, BigMHC, and DeepImmuno). Standalone IEDB is installed at
``/opt/iedb`` (``--iedb-install-directory /opt/iedb``).

An example on how to run pVACseq using Docker can be found on the :ref:`Getting Started <pvacseq_docker>` page.

Common Workflow Language (CWL) and Workflow Description Language (WDL) tool wrappers for pVACseq, pVACfuse, and
pVACvector (CWL only) can be downloaded using the ``pvactools download_cwls`` and
``pvactools download_wdls`` commands, repectively.

Download CWL tool wrappers
__________________________

.. program-output:: pvactools download_cwls -h

Download WDL tool wrappers
__________________________

.. program-output:: pvactools download_wdls -h
