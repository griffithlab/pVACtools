Installation
============

pVAC-Seq requires Python 3.5. Before running any installation steps check the Python version installed on your system:

.. code-block:: none

   python -V

If you don't have Python 3.5 installed, we recommend using `Conda <http://conda.pydata.org/docs/py2or3.html>`_ to emulate a Python 3.5. environment. We've encountered problems with users that already have Python 2.* installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.* installed we **strongly** recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use ``pip`` to install pVAC-Seq. Make sure you have ``pip`` installed. ``pip`` is generally included in python distributions, but may need to be upgraded before use. See the `instructions <https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel>`_ for installing or upgrading pip.

After you have pip installed, type the following command on your Terminal(for Mac and Linux users) or the command prompt (for Windows users):

.. code-block:: none

   pip install pvacseq

You can check that pvacseq has been installed under the default environment by listing all installed packages:

.. code-block:: none

   pip list

pip will fetch and install pVAC-Seq and its dependencies for you. After installing, you can run pvacseq directly from the Terminal/command prompt.

If you have an old version pVAC-Seq installed you might want to consider upgrading to the latest version:

.. code-block:: none

   pip install pvacseq --upgrade
