# pVACtools

pVACtools is a cancer immunotherapy suite consisting of the following tools:

**pVACseq**

A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a list of tumor mutations.

**pVACfuse**

A tool for detecting neoantigens resulting from gene fusions.

**pVACvector**

A tool designed to aid specifically in the construction of DNA vector-based cancer vaccines.

**pVACviz**

A browser-based user interface that assists users in launching, managing, reviewing, and visualizing the results of pVACtools processes.

**pVACapi**

The pVACapi provides a HTTP REST interface to the pVACtools suite.

## Citation
Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.

## Installation
pVACtools requires Python 3.5. Before running any installation steps check the Python version installed on your system:

`python -V`

If you don’t have Python 3.5 installed, we recommend using Conda to emulate a Python 3.5. environment. We’ve encountered problems with users that already have Python 2.x installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.x installed we strongly recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use pip to install pVACtools. Make sure you have pip installed. pip is generally included in python distributions, but may need to be upgraded before use. See the instructions for installing or upgrading pip.

After you have pip installed/upgraded, type the following command on your Terminal:

`pip install pvactools`

You can check that pVACtools has been installed under the default environment by listing all installed packages:

`pip list`

You can also check the installed version:

`pvactools -v`

`pip` will fetch and install pVACtools and its dependencies for you. After installing, each tool of the pVACtools suite is available with its own command line tree directly from the Terminal.

If you have an old version of pVACtools installed you might want to consider upgrading to the latest version:

`pip install pvactools --upgrade`

## Documentation

The pVACtools documentation can be found on <a href="http://pvactools.readthedocs.io/">ReadTheDocs</a>.

## Contact
Bug reports or feature requests can be submitted on the <a href="https://github.com/griffithlab/pVACtools/issues">pVACtools Github page</a>. You may also contact us by email at help@pvactools.org.
