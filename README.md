# pVAC-Seq
Cancer immunotherapy has gained significant momentum from recent clinical successes of checkpoint blockade inhibition. Massively parallel sequence analysis suggests a connection between mutational load and response to this class of therapy. Methods to identify which tumor-specific mutant peptides (neoantigens) can elicit anti-tumor T cell immunity are needed to improve predictions of checkpoint therapy response and to identify targets for vaccines and adoptive T cell therapies. Here, we provide a cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq).
http://www.genomemedicine.com/content/8/1/11

## New in version 3.0.5
<ul>
<li>Bugfix: The generation of the fasta file would fail for some insertions with a range position. This is now fixed.</li>
<li>Bugfix: The generation of the fasta file would fail if the wildtype or downstream sequences were too long. The size limit for these fields has been increased to the user system's maximum supported size. This error might still occur if the sequences are longer than that.</li>
<li>Bugfix: When rerunning a command an error would occur if the <code>tmp</code> subdirectory already exists in the output directory. This has now been fixed.</li>
</ul>

## New in version 3.0.4
<ul>
<li>Certain intermediate files are now written into a <code>tmp</code> directory underneath the main output directory. This <code>tmp</code> directory will be deleted at the end of a successful run unless the <code>--keep-tmp-files</code> flag is set.</li>
<li>Intermediate files will now not be reprocessed if they already exist in the output directory. This can be helpful if a run exits early, for example, when a 500 Error was returned by IEDB. In this case the user can now simply run the same <code>pvacseq run</code> command again and the run will pick up where it failed previously.</li>
<li>We added a new option <code>--fasta-size</code> that the user can set to specify how many FASTA entries at a time will be included in a request to the IEDB RESTful API. The default is 200 but certain variants or prediction algorithms might warrant a smaller number of FASTA entries in order to avoid timeouts from IEDB.</li>
<li>Bugfix: The generation of the fasta file would fail for frameshift mutations with a range position. This is now fixed.</li>
<li>Bugfix: Previously a run might fail if certain intermediate files weren't created.</li>
<li>Bugfix: Using <code>.</code> in the output directory name and the sample name would previously result in errors. This has now been fixed.</li>
<li>Bugfix: Using a relative directory path for the output directory would previsouly result in an error. This is now fixed.</li>
</ul>

## New in version 3.0.3
<ul>
<li>Bugfix: The binding filter used to filter out all but the top peptide candidate for a variant even if the <code>--top-result-per-mutation</code> flag wasn't set. This is now fixed and the top-result-per-mutation filtering only happens when the flag is set.</li>
<li>Bugfix: For large input files the mutant protein sequence wasn't being correctly matched to a wildtype protein sequence. This issue has been corrected.</li>
</ul>

## New in version 3.0.2
<ul>
<li>Bugfix: Some allele names in the list of valid alleles were incorrect. The list has been updated.</li>
<li>If the generate_fasta step creates an empty file during the execution of a run the run will terminate early.</li>
</ul>

## New in version 3.0.0
<ul>
<li>pVAC-Seq now uses the IEDB RESTful interface for making epitope binding predictions. A local install of NetMHC3.4 is no longer required. By using IEDB the user now has a choice between several prediction algorithms, including NetMHC (3.4), NetMHCcons (1.1), NetMHCpan (2.8), PickPocket (1.1), SMM, and SMMPMBEC.</li>
<li>The user can now set the <code>--top-result-per-mutation</code> flag in order to only output the top scoring candidate per allele-length per mutation.</li>
<li>Since it is now possible to run multiple epitope prediction algorithms at the same time, the scores for each candidate epitope are aggregated as <code>Median MT score All Methods</code>, which is the median mutant ic50 binding score of all chosen prediction methods, and the <code>Best MT score</code>, which is the lowest mutant ic50 binding score of all chosen prediction methods. For the Best MT score we also output the <code>Corresponding WT score</code> and the <code>Best MT Score Method</code>. Individual ic50 binding score for each prediction method are also outputted. The user can specify which metric to use for filtering by setting the <code>--top-score-metric</code> argument to either <code>lowest</code> or <code>median</code>.</li>
</ul>

## New in version 2.0.2
<ul>
<li>Bugfix: There was a problem in version 2.0.1 where pVAC-Seq would hang while calling NetMHC under certain cirumstances. This is now fixed.</li>
<li>Bugfix: When using multiple alleles or epitope lengths, pVAC-Seq would not output all candidate epitopes after running the binding filter. This has now been fixed.</li>
</ul>

## New in version 2.0.1
<ul>
<li>Silence output from NetMHC</li>
</ul>

## New in version 2.0.0
<ul>
<li>Supports inframe indels and frameshifts.</li>
<li>Supports VCF as the input file format.</li>
</ul>

## Citation
Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.

## Installation Instructions
pVAC-Seq requires Python 3.5. You can check your installed Python version by running `python -V`. If you don't have Python 3.5 installed, we recommend using <a href="http://conda.pydata.org/docs/py2or3.html">Conda</a> to emulate a Python 3.5. environment. We've encountered problems with users that already have Python 2.x installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.x installed we <b>strongly</b> recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use `pip` to install pVAC-Seq. Make sure you have `pip` installed.  `pip` is generally included in python distributions, but may need to be upgraded before use.  See the <a href="https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel">instructions</a> for installing or upgrading pip.

After you have pip installed, type the following command on your Terminal(for Mac and Linux users) or the command prompt (for Windows users): `pip install pvacseq`. You can check that pvacseq has been installed under the default environment by running `pip list`.

pip will fetch and install pVAC-Seq and its dependencies for you.  After installing, you can run `pvacseq` directly from the Terminal/command prompt.

If you already have pVAC-Seq installed but would like to upgrade to the latest version, you can do so by running `pip install pvacseq --upgrade`.

##Prerequisites
###<b>VEP</b>

The input to the pVAC-Seq pipeline is a VEP annotated VCF. In addition to the standard VEP annotations, pVAC-Seq also requires the annotations provided by the Downstream and Wildtype VEP plugins. To create a VCF for use with pVAC-Seq follow these steps:
- Download and install the VEP command line tool following the instructions found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/index.html">here</a>.
- Download the VEP_plugins from their <a href="https://github.com/Ensembl/VEP_plugins">Github repository</a>.
- Copy the Wildtype plugin provided with the pVAC-Seq package to the folder with the other VEP_plugins by running `pvacseq install_vep_plugin`.
- Run VEP on the input vcf with at least the following options:<br>
`--format vcf`<br>
`--vcf`<br>
`--symbol`<br>
`--plugin Downstream`<br>
`--plugin Wildtype`<br>
`--terms SO`<br>
The `--dir_plugins <VEP_plugins directory>` option may need to be set depending on where the VEP_plugins were installed to. Additional VEP options that might be desired can be found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html">here</a>.<br>
<b>Example VEP Command</b><br>
`perl variant_effect_predictor.pl --input_file <input VCF> --format vcf --output_file <output VCF> --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype [--dir_plugins <VEP_plugins directory>]`

## Pipeline Overview
![alt text][overview]
[overview]:
https://raw.githubusercontent.com/wiki/griffithlab/pVAC-Seq/images/pvacseq-code-pythonv3.0.0.png

## pvacseq commands
### run
`pvacseq run <input VCF> <sample name> <allele name> <epitope length> <prediction_algorithm> <output directory> [-l peptide sequence length] [--top-result-per-mutation] [-m top score metric] [-b binding threshold] [-c minimum fold change] [-s fasta size] [--keep-tmp-files]`<br>
Use this command to run the full pVAC-Seq pipeline.  This will internally call the other commands, passing data between them to generate an output TSV file of neoepitope predictions. Multiple alleles and epiope length can be specified as comma-separated lists.

<b>Required inputs</b><br>
<ul>
<li><code>input VCF</code>: A VEP-annotated VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information. (Please see above for instructions)</li>
<li><code>sample name</code>: The name of the sample being processed. This will be used as a prefix for output files.</li>
<li><code>allele name</code>: Name of the allele to use for epitope prediction. Mutliple alles can be specified using a comma-separated list.</li>
<li><code>epitope length</code>: This refers to the length of subpeptides (neoepitopes) to predict. The pipeline can handle multiple lengths, which can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.</li>
<li><code>prediction algorithm</code>: The prediction algorithm to use. The available choices are <code>NetMHC</code>, <code>NetMHCcons</code>, <code>NetMHCpan</code>, <code>PickPocket</code>, <code>SMM</code>, and <code>SMMPMBEC</code>. Multiple prediction algorithms can be specified, separated by spaces.</li>
<li><code>Output directory</code>: The directory for writing all result files.</li>
</ul>

<b>Optional inputs</b><br>
<ul>
<li><code>peptide sequence length</code>: Length of the peptide sequence to use when creating the FASTA. See "Additional Information" for details. This is set to 21 by default.
<li><code>top result per mutation</code>: When this flag is set only the top scoring candidate per allele-length per mutation will be outputted. By default this is set to false)</li>
<li><code>top score metric</code>: The user can chose which ic50 scoring metric to will be used when filtering epitopes
<ul>
<li>lowest: Best MT Score - lowest mutant ic50 binding score of all chosen prediction methods.</li>
<li>median: Median MT Score All Methods - median mutant ic50 binding score of all chosen prediction methods.</li>
</ul>
By default this argument is set to median.</li>
<li><code>binding threshold</code>: The user can choose to report only epitopes where the mutant allele has IC50 binding scores below this value. By default, pVAC-Seq uses a cutoff of 500.
<li><code>minimum fold change</code>: This parameter sets the minimum fold change between mutant binding score and wild-type score to use for filtering. The default is 0, which filters no results. Using 1 will require that binding is better to the MT than WT.</li>
<li><code>fasta size</code>: The user can specify the number of fasta entries that will be submitted to the IEDB RESTful API at a time. The default is 200 but certain variants or prediction algorithms might warrant a smaller number of FASTA entries in order to avoid timeouts from IEDB.</li>
<li><code>--keep-temp-files</code>: When this flag is set the <code>tmp</code> directory and the intermediate files it contains will not get deleted after a succcessful run. This might be useful for debugging purposes.</li>
</ul>

### convert_vcf
`pvacseq convert_vcf <input VCF> <output TSV file>`<br>
Run this command to generate a TSV file with annotated variants from a VEP-annotated VCF.

### generate_fasta
`pvacseq generate_fasta <input TSV file> <peptide sequence length> <output FASTA file>`<br>
Run this command to generate a FASTA file for wildtype(WT) and mutant(MT) amino acid sequences for MHC Class I epitope prediction. The length of the amino acid sequences is determined by the peptide sequence length specified. The input file is the properly formatted TSV file of annotated variants.

### generate_fasta_key
`pvacseq generate_fasta_key <input FASTA file> <output key file>`<br>
IEDB strips off the name of the FASTA header. This command generates a key file to lookup each IEDB output entry to its original entry in the FASTA file.

## call_iedb
`pvacseq call_iedb <input FASTA file> <output IEDB file> <IEDB analysis method> <allele> <epitope length>`<br>
This command make epitope binding predicitions using the IEDB RESTful interface and writes the result to a file.

### parse_output
`pvacseq parse_output <IEDB files> <input TSV file> <input key file> <output parsed TSV file> [--top-result-per-mutation] [-m <lowest|median>]`<br>
After running IEDB, this command parses the output from the IEDB RESTful API calls. It combines the IEDB output files for multiple prediction algorithms that have the same allele and epitope lengths. It uses a special key file to link each IEDB result entry to the original entry from the input TSV file. The parsed TSV output file contains predictions for the mutant as well as the wildtype version of the epitope, and compares binding affinities for the same. When multiple prediction algorithms are used the parser will find the best mutant ic50 score as well as the median mutant ic50 score. The file also contains gene and transcript information from the input TSV file.

### combine_parsed_outputs
`pvacseq combine_parsed_outputs <input parsed TSV file> <output combined parsed TSV file>`<br>
Combines all parsed output IEDB files into one file. Each parsed output IEDB file contains entries for the same allele and epitope length. This step combines parsed files from multiple alleles and epitope lengths into one single output TSV file.

### binding_filter
`pvacseq binding_filter <input combined parsed TSV file> <output filtered TSV file> [-b binding threshold] [-c minimum fold change] [-m <lowest|median>]`<br>
Takes combined parsed epitope file for different allele-length combinations and outputs best candidates per gene based on binding affinities.

### coverage_filter
`pvacseq coverage_filters <input TSV file> <output filtered TSV file> [--normal-cov normal coverage cutoff] [--tdna-cov tumor DNA coverage cutoff] [--trna-cov tumor RNA coverage cutoff] [--normal-vaf normal vaf cutoff] [--tdna-vaf tumor DNA vaf cutoff] [--trna-vaf tumor RNA vaf cutoff] [--expn-val gene expression (fpkm) cutoff]`<br>
Depending on the type(s) of sequencing data available, a variety of coverage and expression based filters can be installed. The input file should contain the predicted epitopes along with read counts appended as additional columns. If specific type of sequencing data is not available, the columns can be left off. Column order is not important but the names of the headers for the columns containing coverage information is. The headers need to be named as follows:<br>
Normal Ref Count<br>
Normal Var Count<br>
Tumor DNA Ref Count<br>
Tumor DNA Var Count<br>
Tumor RNA Ref Count<br>
Tumor RNA Var Count<br>
Gene Exp FPKM<br>

### download_example_data
`pvacseq download_example_data <destination directory>`
Downloads a set of example data files to the directory specififed.

### install_vep_plugin
`pvacseq install_vep_plugin <vep plugins path>`
Installs the Wildtype VEP plugin into the specified directory.

### valid_alleles
`pvacseq valid_alleles [-p <prediction_algorithm>]`<br>
Shows a list of valid allele names. If the `-p` option is specified with a prediction algorithm than only the alleles available for that predicion algorithm will be displayed. `prediction_algorithm` can be one of `NetMHC`, `NetMHCcons`, `NetMHCpan`, `PickPocket`, `SMM`, or `SMMPMBEC`.

## Additional Information

Since the goal of the pVAC-Seq pipeline is to predict putative 'neo'antigens, we only consider a sub-section of the transcript sequence encompassing the mutated amino acid.
<ol type="a">
<li>In the following figure, the amino acid FASTA sequence is built using 10 flanking amino acids on each side of the mutated amino acid. The preceding or succeeding 20 amino acids are taken if the mutation lies near the end or beginning of the transcript, respectively.</li>
<li>All predicted candidate peptides from epitope prediction software based on selected k-mer window size.</li>
<li>Only localized peptides (those containing the mutant amino acid) are considered to compare to wild-type counterpart.</li>
</ol>
 ![alt text][logo]
 [logo]:
 https://github.com/jhundal/src/blob/master/bin/images/Fig1_fastav2.png
