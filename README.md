# pVAC-Seq 
A cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq) 

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.
## Pipeline Overview
![alt text][overview]
[overview]:
https://github.com/jhundal/src/blob/master/bin/images/pvacseq-code.jpg

## Installation Instructions
pVAC-Seq could either be used via cloning the git repo or by downloading the tarball.
We strongly recommend the users to access it via git repo since updates are pushed before they make it to the stable release.

 1. <b> OPTION 1 : Clone pVAC-Seq git repository: </b>
        Type the following command on your Terminal(for Mac and Linux users) or the command prompt (for Windows users):
   `git clone git@github.com:griffithlab/pVAC-Seq.git  ` 
  
 For more information, follow instructions for cloning a git repo :                               				  https://help.github.com/articles/cloning-a-repository/

 2. <b> OPTION 2 : Download and decompress the tar ball</b>: 

 * `gunzip pVAC-Seq-1.0.0-beta.tar.gz`
 * `tar -xzf pVAC-Seq-1.0.0-beta.tar.gz`
 * Make sure you see the following directories:
 
                        LICENSE
                        README.md
                        bin/
                        pVAC-Seq.pl
                        test_data/`

 * Adjust the perl shebang line of each .pl and .sh script in the bin/ folder as needed
 * Execute(run) pVAC-Seq pipeline by typing the following command and providing necessary inputs :
 
 	`./pVAC-Seq.pl`


##Pre-requisites
###<b>NetMHC </b>

Since we use NetMHC to predict binding affinities, it is one of the major pre-requisites to run pVAC-Seq.pl.
Please download NetMHC using appropriate licensing. Download and installation instructions are provided on the website:
  http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHC

   Also note that there might be a minor bug in the NetMHC package, which can be fixed as follows :
	 http://www.cbs.dtu.dk/services/NetMHC/correct_bug.php
	 
  Once NetMHC is properly installed and tested, pVAC-Seq.pl expects the path to the installtion directory.


## Inputs for pVAC-Seq.pl

1. <b>NetMHC installation path: </b> Provide path to the NetMHC installation directory (please see above for installation instructions)

2. <b>TSV file of annotated variants</b>: The program expects an annotated file of variants in tab-separated format. Any choice of aligner or variant caller can be used. The following columns are expected as part of the TSV file (in the same order) along with the header row:

 |chromosome_name| start | stop | reference | variant | gene_name | transcript_name | amino_acid_change | ensembl_gene_id |wildtype_amino_acid_sequence
--- | --- | --- | ---| ---| ---| ---| ---| ---| ---| ---| ---| ---
1	| 92163648|	92163648|	G	|A	|TGFBR3	|ENST00000212355	|P776L	|ENSG00000069702|	MTSHYVIAIFALMSSCLATAGPEPGALCELSPVSASHPVQALMESFTVLSGCASRGTTGLPQEVHVLNLRTAGQGPGQLQREVTLHLNPISSVHIHHKSVVFLLNSPHPLVWHLKTERLATGVSRLFLVSEGSVVQFSSANFSLTAETEERNFPHGNEHLLNWARKEYGAVTSFTELKIARNIYIKVGEDQVFPPKCNIGKNFLSLNYLAEYLQPKAAEGCVMSSQPQNEEVHIIELITPNSNPYSAFQVDITIDIRPSQEDLEVVKNLILILKCKKSVNWVIKSFDVKGSLKIIAPNSIGFGKESERSMTMTKSIRDDIPSTQGNLVKWALDNGYSPITSYTMAPVANRFHLRLENNAEEMGDEEVHTIPPELRILLDPGALPALQNPPIRGGEGQNGGLPFPFPDISRRVWNEEGEDGLPRPKDPVIPSIQLFPGLREPEEVQGSVDIALSVKCDNEKMIVAVEKDSFQASGYSGMDVTLLDPTCKAKMNGTHFVLESPLNGCGTRPRWSALDGVVYYNSIVIQVPALGDSSGWPDGYEDLESGDNGFPGDMDEGDASLFTRPEIVVFNCSLQQVRNPSSFQEQPHGNITFNMELYNTDLFLVPSQGVFSVPENGHVYVEVSVTKAEQELGFAIQTCFISPYSNPDRMSHYTIIENICPKDESVKFYSPKRVHFPIPQADMDKKRFSFVFKPVFNTSLLFLQCELTLCTKMEKHPQKLPKCVPPDEACTSLDASIIWAMMQNKKTFTKPLAVIHHEAESKEKGPSMKEPNPISPPIFHGLDTLTVMGIAFAAFVIGALLTGALWYIYSHTGETAGRQQVPTSPPASENSSAAHSIGSTQSTPCSSSSTA|
1	|108291655|	108291655|	C|	T|	VAV3|	ENST00000490388|	G474D|	ENSG00000134215|	XCAQWLIHCKVLPTNHRVTWDSAQVFDLAQTLRDGVLLCQLLNNLRAHSINLKEINLRPQMSQFLCLKNIRTFLTACCETFGMRKSELFEAFDLFDVRDFGKVIETLSRLSRTPIALATGIRPFPTEESINDEDIYKGLPDLIDETLVEDEEDLYDCVYGEDEGGEVYEDLMKAEEAHQPKCPENDIRSCCLAEIKQTEEKYTETLESIEKYFMAPLKRFLTAAEFDSVFINIPELVKLHRNLMQEIHDSIVNKNDQNLYQVFINYKERLVIYGQYCSGVESAISSLDYISKTKEDVKLKLEECSKRANNGKFTLRDLLVVPMQRVLKYHLLLQELVKHTTDPTEKANLKLALDAMKDLAQYVNEVKRDNETLREIKQFQLSIENLNQPVLLFGRPQGDGEIRITTLDKHTKQERHIFLFDLAVIVCKRKGDNYEMKEIIDLQQYKIANNPTTDKENKKWSYGFYLIHTQGQNGLEFYCKTKDLKKKWLEQFEMALSNIRPDYADSNFHDFKMHTFTRVTSCKVCQMLLRGTFYQGYLCFKCGARAHKECLGRVDNCGRVNSGEQGTLKLPEKRTNGLRRTPKQVDPDVPCLLHFFISMAPATRSIVKSQKKNKKF

  Any annotation database could be used for providing this information, as long as gene id, transcript id and wildtype transcript sequence is provided. 

3. <b>Variant peptide sequence length</b>: Since the goal of the pVAC-Seq pipeline to predict putative 'neo'antigens, we only consider a sub-section of the transcript sequence encompassing the mutated amino acid. 

     (a) In the following figure, amino acid FASTA sequence is built using 10 flanking amino acids on each side of the mutated amino acid. The preceding or succeeding 20 amino acids are taken if the mutation lies near the end or beginning of the transcript, respectively.

     (b). All predicted candidate peptides from epitope prediction software based on selected k-mer window size. 

     (c). Only localized peptides (those containing the mutant amino acid) are considered to compare to wild-type counterpart. 

     (d). The ‘best candidate’ (lowest MT binding score) per mutation is chosen across all specified k-mers that were used as input.
![alt text][logo]
[logo]:
https://github.com/jhundal/src/blob/master/bin/images/Fig1_fastav2.png
4.  <b>Epitope length</b> : This refers to the length of subpeptides(neoepitopes) to predict. The pipeline can handle multiple lengths that can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.

5. <b> Binding-cutoff </b> : The user can choose to report only epitopes where the mutant allele has IC50 binding scores below this value. By default, we recommend choosing high to medium binding epitopes and use a cutoff of 500.
6. <b> Minimum Fold Change (min-fc):</b> This parameter is used to set the minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT). 

## Individual Modules in bin/

1. <b> GenerateVariantSequences.pl </b>: Run this script to generate a FASTA file for wildtype(WT) and mutant(MT) 21-mer amino acid sequences for MHC Class I epitope prediction. The input file is the properly formatted TSV file of annotated variants. The following columns are expected as part of the TSV file (in the same order) along with the header row:

 |chromosome_name| start | stop | reference | variant | gene_name | transcript_name | amino_acid_change | ensembl_gene_id |wildtype_amino_acid_sequence
--- | --- | --- | ---| ---| ---| ---| ---| ---| ---| ---| ---| ---
1	| 92163648|	92163648|	G	|A	|TGFBR3	|ENST00000212355	|P776L	|ENSG00000069702|	MTSHYVIAIFALMSSCLATAGPEPGALCELSPVSASHPVQALMESFTVLSGCASRGTTGLPQEVHVLNLRTAGQGPGQLQREVTLHLNPISSVHIHHKSVVFLLNSPHPLVWHLKTERLATGVSRLFLVSEGSVVQFSSANFSLTAETEERNFPHGNEHLLNWARKEYGAVTSFTELKIARNIYIKVGEDQVFPPKCNIGKNFLSLNYLAEYLQPKAAEGCVMSSQPQNEEVHIIELITPNSNPYSAFQVDITIDIRPSQEDLEVVKNLILILKCKKSVNWVIKSFDVKGSLKIIAPNSIGFGKESERSMTMTKSIRDDIPSTQGNLVKWALDNGYSPITSYTMAPVANRFHLRLENNAEEMGDEEVHTIPPELRILLDPGALPALQNPPIRGGEGQNGGLPFPFPDISRRVWNEEGEDGLPRPKDPVIPSIQLFPGLREPEEVQGSVDIALSVKCDNEKMIVAVEKDSFQASGYSGMDVTLLDPTCKAKMNGTHFVLESPLNGCGTRPRWSALDGVVYYNSIVIQVPALGDSSGWPDGYEDLESGDNGFPGDMDEGDASLFTRPEIVVFNCSLQQVRNPSSFQEQPHGNITFNMELYNTDLFLVPSQGVFSVPENGHVYVEVSVTKAEQELGFAIQTCFISPYSNPDRMSHYTIIENICPKDESVKFYSPKRVHFPIPQADMDKKRFSFVFKPVFNTSLLFLQCELTLCTKMEKHPQKLPKCVPPDEACTSLDASIIWAMMQNKKTFTKPLAVIHHEAESKEKGPSMKEPNPISPPIFHGLDTLTVMGIAFAAFVIGALLTGALWYIYSHTGETAGRQQVPTSPPASENSSAAHSIGSTQSTPCSSSSTA|
1	|108291655|	108291655|	C|	T|	VAV3|	ENST00000490388|	G474D|	ENSG00000134215|	XCAQWLIHCKVLPTNHRVTWDSAQVFDLAQTLRDGVLLCQLLNNLRAHSINLKEINLRPQMSQFLCLKNIRTFLTACCETFGMRKSELFEAFDLFDVRDFGKVIETLSRLSRTPIALATGIRPFPTEESINDEDIYKGLPDLIDETLVEDEEDLYDCVYGEDEGGEVYEDLMKAEEAHQPKCPENDIRSCCLAEIKQTEEKYTETLESIEKYFMAPLKRFLTAAEFDSVFINIPELVKLHRNLMQEIHDSIVNKNDQNLYQVFINYKERLVIYGQYCSGVESAISSLDYISKTKEDVKLKLEECSKRANNGKFTLRDLLVVPMQRVLKYHLLLQELVKHTTDPTEKANLKLALDAMKDLAQYVNEVKRDNETLREIKQFQLSIENLNQPVLLFGRPQGDGEIRITTLDKHTKQERHIFLFDLAVIVCKRKGDNYEMKEIIDLQQYKIANNPTTDKENKKWSYGFYLIHTQGQNGLEFYCKTKDLKKKWLEQFEMALSNIRPDYADSNFHDFKMHTFTRVTSCKVCQMLLRGTFYQGYLCFKCGARAHKECLGRVDNCGRVNSGEQGTLKLPEKRTNGLRRTPKQVDPDVPCLLHFFISMAPATRSIVKSQKKNKKF

  Any annotation database could be used for providing this information, as long as gene id, transcript id and wildtype transcript sequence is provided. 
2. <b>GenerateFastaKey.pl</b>: NetMHC strips off the name of the FASTA header that contains gene names and type of sequence (WT vs MT). This module generates a key file to lookup original gene names in the output file of NetMHC 3.4 from the original 21-mer FASTA file for wildtype(WT) and mutant(MT) proteins.

3. <b>ParseOutputNetmhc.pl</b>:  After running NETMHC3.4, this module parses the output for MHC Class I epitope prediction. It uses a special key file generated that could be generated using GenerateFastaKey.pl.The parsed TSV file contains predictions for the mutant as well as the wildtype version of the epitope, and compares binding affinities for the same. 

4. <b>BindingFilter.pl</b>: Takes in a file of files with path to parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinities.

5. <b>CoverageFilters.pl</b>: Depending on the type(s) of sequencing data available, a variety of coverage and expression based filters could be used. The input file should contain the predicted epitopes along with read counts appended as additional columns. <b>Please note that if specific type of sequencing data is not available, the user should enter n/a in those columns, and set appropriate flags while running the script.</b> Column order should be preserved.

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

6. <b>GenerateFastaForNetChop.pl</b>: Takes in the filtered file generated from pVAC-Seq.pl, BindingFilter.pl or CoverageFilters.pl, and generates a FASTA file of MT epitope sequences. This FASTA file can be used as an input to NetChop to evaluate predictions for cleavage sites of the human proteasome.
