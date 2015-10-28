####Example command : 

./pVAC-Seq.pl  --input-file ~/pVAC-Seq/test_data/annotated_variants.tsv --sample-name Test --netmhc-path ~/netMHC-3.4/netMHC --allele HLA-A29:02 --varpeptide-length 21 --epitope-length 9 --output-dir ~/pVAC-Seq/test_data/


####Command Explanation:
`./pVAC-Seq.pl`
`--input-file ~/pVAC-Seq/test_data/annotated_variants.tsv`
#######Use the TSV as an example to format your input file
`--sample-name Test` 
#######Specify your preferred sample name ; this will used as a prefix for all output files
`--netmhc-path ~/netMHC-3.4/netMHC `
####### Install NetMHC per installation instruction and provide path to the directory
`--allele HLA-A29:02 `
####### One allele, or comma-separated list of different alleles
`--varpeptide-length 21 `
####### Length of FASTA sequence with mutated amino acid
`--epitope-length 9 `
####### Predicted epitope length ; Comma-separated list may be used
`--output-dir ~/pVAC-Seq/test_data/`
####### Directory where all result files will be written.
