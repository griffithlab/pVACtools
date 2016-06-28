# pVAC-Seq
A cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq)
http://www.genomemedicine.com/content/8/1/11

## Citation
Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.
## Pipeline Overview
![alt text][overview]
[overview]:
https://github.com/jhundal/src/blob/master/bin/images/pvacseq-code.jpg

## Installation Instructions
  Make sure you have pip installed.  pip is generally included in python distributions, but may need to be upgraded before use.  See these instructions for installing or upgrading pip: https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel

  After you have pip installed, type the following command on your Terminal(for Mac and Linux users) or the command prompt (for Windows users): `pip install pvacseq`

  pip will fetch and install pVAC-Seq and its dependencies for you.  After installing, you can run `pvacseq` directly from the Terminal/command prompt

##Prerequisites
###<b>VEP</b>

The input to the pVAC-Seq pipeline is a VEP annotated VCF. In addition to the standard VEP annotations, pVAC-Seq also requires the annotations provided by the Downstream and Wildtype VEP plugins. To create a VCF for use with pVAC-Seq follow these steps:
- Download and install the VEP command line tool following the instructions found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/index.html">here</a>.
- Download the VEP_plugins from their <a href="https://github.com/Ensembl/VEP_plugins">Github repository</a>
- Copy the Wildtype plugin provided with the pVAC-Seq package to the folder with the other VEP_plugins by running `pvacseq install_vep_plugin`.
- Run VEP on the input vcf with at least the following options:<br>
`--format vcf`<br>
`--vcf`<br>
`--symbol`<br>
`--plugin Downstream`<br>
`--plugin Wildtype`<br>
`--terms SO`<br>
The `--dir_plugins <VEP_plugins directory>` option may need to be set depending on where the VEP_plugins were installed to. Additional VEP options that might be desired can be found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html"here</a>.

###<b>NetMHC 3.4</b>

pVAC-Seq uses NetMHC 3.4 to predict binding affinities. NetMHC 3.4 can be downloaded <a href="http://www.cbs.dtu.dk/cgi-bin/sw_request?netMHC+3.4">here</a>. Once NetMHC is properly installed and tested, pVAC-Seq expects the path to the installtion directory.


## pvacseq commands
<b>run</b><br>
`pvacseq run <input VCF file> <sample name> <NetMHC installation path> <allele name> <epitope length> <ouput directory> [-l peptide sequence length] [-b binding threshold] [-c minimum fold change]`<br>
Run this command to automate the pVAC-Seq pipeline.  This will internally call the other commands, passing data between them to generate an TSV file of neoepitope predictions. Multiple alleles and epiope length can be specified as comma-separated lists.

<b>convert_vcf</b><br>
`pvacseq convert_vcf <input VCF file> <output TSV file>`
Run this command to generate a TSV file with annotated variants from a VEP-annotated VCF.

<b>generate_fasta</b><br>
`pvacseq generate_fasta <input TSV file> <peptide sequence length> <output FASTA file>`<br>
Run this command to generate a FASTA file for wildtype(WT) and mutant(MT) amino acid sequences for MHC Class I epitope prediction. The length of the amino acid sequences is determined by the peptide sequence length specified. The input file is the properly formatted TSV file of annotated variants.

<b>generate_fasta_key</b><br>
`pvacseq generate_fasta_key <input FASTA file> <output key file>`<br>
NetMHC strips off the name of the FASTA header. This command generates a key file to lookup each NetMHC output entry to its original entry in the FASTA file.

<b>parse_output</b><br>
`pvacseq parse_output <NetMHC output file> <input TSV file> <FASTA key file> <output parsed file>`<br>
After running NetMHC 3.4, this command parses the output for MHC Class I epitope prediction. It uses a special key file to link each NetMHC result entry to the original entry from the input TSV file. The parsed TSV output file contains predictions for the mutant as well as the wildtype version of the epitope, and compares binding affinities for the same. It also contains gene and transcript information from the input TSV file.

<b>binding_filter</b><br>
`pvacseq binding_filter <input TSV file> <output file> [-b binding threshold] [-c minimum fold change]`<br>
Takes a comma-separated list of parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinities.

<b>coverage_filter</b><br>
`pvacseq coverage_filters <input TSV file> <output file> [--normal-cov normal coverage cutoff] [--tdna-cov tumor DNA coverage cutoff] [--trna-cov turmor DNA coverage cutoff] [--normal-vaf normal vaf cutoff] [--tdna-vaf tumor DNA vaf cutoff] [--trna-vaf tumor RNA vaf cutoff] [--expn-val gene expression (fpkm) cutoff]`<br>
Depending on the type(s) of sequencing data available, a variety of coverage and expression based filters can be installed. The input file should contain the predicted epitopes along with read counts appended as additional columns. If specific type of sequencing data is not available, the columns can be left off. Column order is not important.

The input TSV file contains the following columns in tab-separated format:<br>
Chromosome<br>
Start<br>
Stop<br>
Reference<br>
Variant<br>
Transcript<br>
Ensembl Gene ID<br>
Variant Type<br>
Mutation<br>
Protein Position<br>
Gene Name<br>
HLA Allele<br>
Peptide Length<br>
Sub-peptide Position<br>
MT score<br>
WT score<br>
MT epitope seq<br>
WT epitope seq<br>
Fold Change<br>
Normal Ref Count<br>
Normal Var Count<br>
Tumor DNA Ref Count<br>
Tumor DNA Var Count<br>
Tumor RNA Ref Count<br>
Tumor RNA Var Count<br>
Gene Exp FPKM<br>

<b>download_example_data</b><br>
`pvacseq download_example_data <destination directory>`
Downloads a set of example data files to the directory specififed.

<b>install_vep_plugin</b><br>
`pvacseq install_vep_plugin <vep plugins path>`
Installs the Wildtype VEP plugin into the specified directory.

## Inputs for pvacseq commands

 1. <b>NetMHC installation path: </b> Provide path to the NetMHC installation directory (please see above for installation instructions)

 2. <b>TSV file of annotated variants</b>: The program expects an annotated file of variants in tab-separated format. Any choice of aligner or variant caller can be installed. The following columns are expected as part of the TSV file (in the same order) along with the header row:

  |chromosome_name| start | stop | reference | variant | gene_name | transcript_name | amino_acid_change | ensembl_gene_id |wildtype_amino_acid_sequence
 --- | --- | --- | ---| ---| ---| ---| ---| ---| ---| ---| ---| ---
 1	| 92163648|	92163648|	G	|A	|TGFBR3	|ENST00000212355	|P776L	|ENSG00000069702|	MTSHYVIAIFALMSSCLATAGPEPGALCELSPVSASHPVQALMESFTVLSGCASRGTTGLPQEVHVLNLRTAGQGPGQLQREVTLHLNPISSVHIHHKSVVFLLNSPHPLVWHLKTERLATGVSRLFLVSEGSVVQFSSANFSLTAETEERNFPHGNEHLLNWARKEYGAVTSFTELKIARNIYIKVGEDQVFPPKCNIGKNFLSLNYLAEYLQPKAAEGCVMSSQPQNEEVHIIELITPNSNPYSAFQVDITIDIRPSQEDLEVVKNLILILKCKKSVNWVIKSFDVKGSLKIIAPNSIGFGKESERSMTMTKSIRDDIPSTQGNLVKWALDNGYSPITSYTMAPVANRFHLRLENNAEEMGDEEVHTIPPELRILLDPGALPALQNPPIRGGEGQNGGLPFPFPDISRRVWNEEGEDGLPRPKDPVIPSIQLFPGLREPEEVQGSVDIALSVKCDNEKMIVAVEKDSFQASGYSGMDVTLLDPTCKAKMNGTHFVLESPLNGCGTRPRWSALDGVVYYNSIVIQVPALGDSSGWPDGYEDLESGDNGFPGDMDEGDASLFTRPEIVVFNCSLQQVRNPSSFQEQPHGNITFNMELYNTDLFLVPSQGVFSVPENGHVYVEVSVTKAEQELGFAIQTCFISPYSNPDRMSHYTIIENICPKDESVKFYSPKRVHFPIPQADMDKKRFSFVFKPVFNTSLLFLQCELTLCTKMEKHPQKLPKCVPPDEACTSLDASIIWAMMQNKKTFTKPLAVIHHEAESKEKGPSMKEPNPISPPIFHGLDTLTVMGIAFAAFVIGALLTGALWYIYSHTGETAGRQQVPTSPPASENSSAAHSIGSTQSTPCSSSSTA|
 1	|108291655|	108291655|	C|	T|	VAV3|	ENST00000490388|	G474D|	ENSG00000134215|	XCAQWLIHCKVLPTNHRVTWDSAQVFDLAQTLRDGVLLCQLLNNLRAHSINLKEINLRPQMSQFLCLKNIRTFLTACCETFGMRKSELFEAFDLFDVRDFGKVIETLSRLSRTPIALATGIRPFPTEESINDEDIYKGLPDLIDETLVEDEEDLYDCVYGEDEGGEVYEDLMKAEEAHQPKCPENDIRSCCLAEIKQTEEKYTETLESIEKYFMAPLKRFLTAAEFDSVFINIPELVKLHRNLMQEIHDSIVNKNDQNLYQVFINYKERLVIYGQYCSGVESAISSLDYISKTKEDVKLKLEECSKRANNGKFTLRDLLVVPMQRVLKYHLLLQELVKHTTDPTEKANLKLALDAMKDLAQYVNEVKRDNETLREIKQFQLSIENLNQPVLLFGRPQGDGEIRITTLDKHTKQERHIFLFDLAVIVCKRKGDNYEMKEIIDLQQYKIANNPTTDKENKKWSYGFYLIHTQGQNGLEFYCKTKDLKKKWLEQFEMALSNIRPDYADSNFHDFKMHTFTRVTSCKVCQMLLRGTFYQGYLCFKCGARAHKECLGRVDNCGRVNSGEQGTLKLPEKRTNGLRRTPKQVDPDVPCLLHFFISMAPATRSIVKSQKKNKKF

   Any annotation database can be installed for providing this information, as long as gene id, transcript id and wildtype transcript sequence are provided.

 3. <b>Variant peptide sequence length</b>: Since the goal of the pVAC-Seq pipeline to predict putative 'neo'antigens, we only consider a sub-section of the transcript sequence encompassing the mutated amino acid.

      (a) In the following figure, the amino acid FASTA sequence is built using 10 flanking amino acids on each side of the mutated amino acid. The preceding or succeeding 20 amino acids are taken if the mutation lies near the end or beginning of the transcript, respectively.

      (b). All predicted candidate peptides from epitope prediction software based on selected k-mer window size.

      (c). Only localized peptides (those containing the mutant amino acid) are considered to compare to wild-type counterpart.

      (d). The ‘best candidate’ (lowest MT binding score) per mutation is chosen across all specified k-mers that were installed as input.
 ![alt text][logo]
 [logo]:
 https://github.com/jhundal/src/blob/master/bin/images/Fig1_fastav2.png
 4.  <b>Epitope length</b> : This refers to the length of subpeptides(neoepitopes) to predict. The pipeline can handle multiple lengths that can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.

 5. <b> Binding-cutoff </b> : The user can choose to report only epitopes where the mutant allele has IC50 binding scores below this value. By default, we recommend choosing high to medium binding epitopes and use a cutoff of 500.

 6. <b> Minimum Fold Change (min-fc):</b> This parameter is installed to set the minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT).

 ####Example command :

 `pvacseq run ~/pVAC-Seq/example_data/annotated_variants.tsv Test ~/netMHC-3.4/netMHC HLA-A29:02 9 ~/pVAC-Seq/example_data/ -l 21`

 ####Command Explanation:
 `pVAC-Seq`
 `~/pVAC-Seq/test_data/annotated_variants.tsv`
 #######Use the TSV as an example to format your input file
 `Test`
 #######Specify your preferred sample name ; this will used as a prefix for all output files
 `~/netMHC-3.4/netMHC `
 ####### Install NetMHC and provide path to the directory
 `HLA-A29:02 `
 ####### One allele, or comma-separated list of different alleles
 `9 `
 ####### Predicted epitope length ; Comma-separated list may be used
 `~/pVAC-Seq/test_data/`
 ####### Directory where all result files will be written.
 `-l 21 `
 ####### Length of FASTA sequence with mutated amino acid
