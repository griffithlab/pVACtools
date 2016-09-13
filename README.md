# pVAC-Seq
pVAC-Seq is a cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. It is being used in studies of checkpoint therapy response and to identify targets for cancer vaccines and adoptive T cell therapies. For more general information, see the <a href="http://www.genomemedicine.com/content/8/1/11">manuscript published in Genome Medicine</a>.

## Citation
Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.

## Installation
pVAC-Seq requires Python 3.5. Before running any installation steps check the Python version installed on your system:

`python -V`

If you don’t have Python 3.5 installed, we recommend using <a href="http://conda.pydata.org/docs/py2or3.html">Conda</a> to emulate a Python 3.5. environment. We’ve encountered problems with users that already have Python 2.x installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.x installed we <b>strongly</b> recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use `pip` to install pVAC-Seq. Make sure you have `pip` installed.  `pip` is generally included in python distributions, but may need to be upgraded before use.  See the <a href="https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel">instructions</a> for installing or upgrading pip.

After you have pip installed, type the following command on your Terminal (for Mac and Linux users) or the Command Prompt (for Windows users):

`pip install pvacseq`

You can check that pvacseq has been installed under the default environment by listing all installed packages:

`pip list`

`pip` will fetch and install pVAC-Seq and its dependencies for you. After installing, you can run `pvacseq` directly from the Terminal/Command Prompt.

If you have an old version of pVAC-Seq installed you might want to consider upgrading to the latest version:

`pip install pvacseq --upgrade`

##Documentation
The pVAC-Seq documentation can be found on <a href="http://pvac-seq.readthedocs.io/">ReadTheDocs</a>.
