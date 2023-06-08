![Test Status](https://github.com/griffithlab/pVACtools/actions/workflows/tests.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/griffithlab/pVACtools/badge.svg?branch=master)](https://coveralls.io/github/griffithlab/pVACtools?branch=master)
[![Docs](https://readthedocs.org/projects/pvactools/badge/?version=latest)](http://pvactools.readthedocs.io/en/latest/?badge=latest)
![External APIs Status](https://github.com/griffithlab/pVACtools/actions/workflows/api_status.yml/badge.svg)
<a href="https://pypi.python.org/pypi/pvactools/">
    <img src="https://img.shields.io/pypi/v/pvactools.svg?maxAge=1000" alt="PyPI" />
</a>

# pVACtools

pVACtools is a cancer immunotherapy suite consisting of the following tools:

**pVACseq**

A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a VCF file.

**pVACbind**

A cancer immunotherapy pipeline for identifying and prioritizing neoantigens from a FASTA file.

**pVACfuse**

A tool for detecting neoantigens resulting from gene fusions.

**pVACvector**

A tool designed to aid specifically in the construction of DNA vector-based cancer vaccines.

**pVACview**

An application based on R Shiny that assists users in reviewing, exploring and prioritizing neoantigens from the results of pVACtools processes for personalized cancer vaccine design.

## Citations
Jasreet Hundal , Susanna Kiwala , Joshua McMichael, Chris Miller, Huiming Xia, Alex Wollam, Conner Liu, Sidi Zhao, Yang-Yang Feng, Aaron Graubert, Amber Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William Gillanders, Elaine Mardis, Obi Griffith, Malachi Griffith. pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens. Cancer Immunology Research. 2020 Mar;8(3):409-420. doi: 10.1158/2326-6066.CIR-19-0401. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/31907209">31907209</a>.

Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, William C. Chapman, Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="https://doi.org/10.1038/s41588-018-0283-9">Accounting for proximal variants improves neoantigen prediction</a>. Nature Genetics. 2018, DOI: 10.1038/s41588-018-0283-9. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/30510237">30510237</a>.

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="https://spdx.org/licenses/BSD-3-Clause-Clear.html">BSD 3-Clause Clear License</a>.

## Installation
pVACtools is written for Linux but some users have been able to run it successfully on Mac OS X. If you are using Windows you will need to set up a Linux environment, for example by setting up a virtual machine.

pVACtools requires Python 3.6 or above. Before running any installation steps, check the Python version installed on your system:

`python -V`

If you don’t have Python 3 installed, we recommend using [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to emulate a Python 3 environment. We’ve encountered problems with users that already have Python 2.x installed when they also try to install Python 3. The defaults will not be set correctly in that case. If you already have Python 2.x installed we **strongly** recommmend using Conda instead of installing Python 3 locally.

Once you have set up your Python 3 environment correctly you can use `pip` to install pVACtools. Make sure you have `pip` installed. `pip` is generally included in python distributions, but may need to be upgraded before use. See the [instructions](https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel) for installing or upgrading `pip`.

After you have `pip` installed, type the following command on your Terminal:

`pip install pvactools`

You can check that `pvactools` has been installed under the default environment like so:

`pip show pvactools`

`pip` will fetch and install pVACtools and its dependencies for you. After installing, each tool of the pVACtools suite is available in its own command line tree directly from the Terminal.

If you have an old version of pVACtools installed you might want to consider upgrading to the latest version:

`pip install pvactools --upgrade`

## Documentation

The pVACtools documentation can be found on <a href="http://pvactools.readthedocs.io/">ReadTheDocs</a>.

## Contact

Bug reports or feature requests can be submitted on the <a href="https://github.com/griffithlab/pVACtools/issues">pVACtools Github page</a>. You may also contact us by email at help@pvactools.org.

## Container images

pVACtools is available as a Docker Image at <a href="https://hub.docker.com/r/griffithlab/pvactools/">DockerHub griffithlab/pvactools</a>.

## Stable release with DOI

[![DOI](https://zenodo.org/badge/102625109.svg)](https://zenodo.org/badge/latestdoi/102625109)
