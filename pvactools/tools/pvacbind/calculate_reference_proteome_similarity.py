import sys

from pvactools.lib.calculate_reference_proteome_similarity import CalculateReferenceProteomeSimilarity

def define_parser():
    return CalculateReferenceProteomeSimilarity.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    CalculateReferenceProteomeSimilarity(
        args.input_file,
        args.input_fasta,
        args.output_file,
        args.match_length,
        args.species,
        'pVACbind',
        args.blastp_path,
        args.blastp_db,
        args.n_threads
    ).execute()

if __name__ == "__main__":
    main()
