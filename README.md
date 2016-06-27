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
###<b>NetMHC </b>

Since we use NetMHC to predict binding affinities, it is one of the major prerequisites to run pVAC-Seq
Once NetMHC is properly installed and tested, pVAC-Seq expects the path to the installtion directory.


## pvacseq commands

1. <b> run </b>:
`pvacseq run <input TSV file> <sample name> <NetMHC installation path> <allele length> <epitope length> <ouput directory> [-l Variant_peptide_sequence_length] [-b Binding_threshold] [-c Minimum_fold_change]`
Run this command to automate the pVAC-Seq pipeline.  This will internally call the other commands, passing data between them to generate an .xls spreadsheet of neoepitope predictions

2. <b> generate_fasta </b>:
`pvacseq generate_fasta <input TSV file> <variant peptide sequence length> <output FASTA file>`
Run this command to generate a FASTA file for wildtype(WT) and mutant(MT) 21-mer amino acid sequences for MHC Class I epitope prediction. The input file is the properly formatted TSV file of annotated variants.
3. <b> generate_fasta_key </b>:
`pvacseq generate_fasta_key <input FASTA file (from generate_fasta)> <output key file>`
NetMHC strips off the name of the FASTA header that contains gene names and type of sequence (WT vs MT). This command generates a key file to lookup original gene names in the output file of NetMHC 3.4 from the original 21-mer FASTA file for wildtype(WT) and mutant(MT) proteins.

4. <b>parse_output</b>:  
`pvacseq parse_output <NetMHC output file> <FASTA key file (from generate_fasta_key)> <output parsed file>`
After running NETMHC3.4, this command parses the output for MHC Class I epitope prediction. It uses a special key file generated that can be generated using generate_fasta_key. The parsed TSV file contains predictions for the mutant as well as the wildtype version of the epitope, and compares binding affinities for the same.

5. <b>binding_filter</b>:
`pvacseq binding_filter <input TSV file> <FOF file containing filepaths to parsed NetMHC files (from parse_output)> <output file> [-b Binding_threshold] [-c Minimum_fold_change]`
Takes in a file of files with path to parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinities.

6. <b>coverage_filters</b>:
`pvacseq coverage_filters *args*`
Depending on the type(s) of sequencing data available, a variety of coverage and expression based filters can be installed. The input file should contain the predicted epitopes along with read counts appended as additional columns. <b>Please note that if specific type of sequencing data is not available, the user should enter n/a in those columns, and set appropriate flags while running the script.</b> Column order should be preserved.

   The Input file contains the following columns in tab-separated format :
   1.	chromosome_name
   2.	start
   3.	stop
   4.	reference
   5.	variant
   6.	gene_name
   7.	transcript_name
   8.	amino_acid_change
   9.	ensembl_gene_id
   10. wildtype_amino_acid_sequence
   11.	GeneName
   12.	HLAallele
   13.	PeptideLength
   14.	SubPeptidePosition
   15.	MTScore
   16.	WTScore
   17.	MTEpitopeSeq
   18.	WTEpitopeSeq
   19.	FoldChange
   20.	NormalRefCount
   21. NormalVarCount
   22. TumorDNARefCount
   23. TumorDNAVarCount
   24. TumorRNARefCount
   25. TumorRNAVarCount
   26. GeneExpFPKM


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
