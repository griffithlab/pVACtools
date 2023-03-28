import sys

from pvactools.lib.calculate_reference_proteome_similarity import CalculateReferenceProteomeSimilarity

def define_parser():
    return CalculateReferenceProteomeSimilarity.parser('pvacfuse')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    CalculateReferenceProteomeSimilarity(
        args.input_file,
        args.input_fasta,
        args.output_file,
        match_length=args.match_length,
        species=args.species,
        file_type='pVACfuse',
        blastp_path=args.blastp_path,
        blastp_db=args.blastp_db,
        peptide_fasta=args.peptide_fasta,
        n_threads=args.n_threads
    ).execute()

if __name__ == "__main__":
    main()
