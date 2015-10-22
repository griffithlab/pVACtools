# pVAC-Seq 
A cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq) 

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.

## Inputs :

1. <b>NetMHC installation path:</b>

  Since we use NetMHC to predict binding affinities, it is one of the major pre-requisites to run pVAC-Seq.pl.
  Please download NetMHC using appropriate licensing. Download and installation instructions are provided on the website:
http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHC

  Also note that there might be a minor bug in the NetMHC package, which can be fixed as follows :
http://www.cbs.dtu.dk/services/NetMHC/correct_bug.php
  
  Once NetMHC is properly installed and tested, pVAC-Seq.pl expects the path to the installtion directory.

2. <b>TSV file of annotated variants</b>:
  The program expects an annotated file of variants in tab-separated format. Any choice of aligner or variant caller can be used. The following columns are expected as part of the TSV file (in the same order) along with the header row:

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
