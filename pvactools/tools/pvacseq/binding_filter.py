import sys

from pvactools.lib.binding_filter import BindingFilter

def define_parser():
    return BindingFilter.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    BindingFilter(
        args.input_file, args.output_file,
        binding_threshold=args.binding_threshold,
        minimum_fold_change=args.minimum_fold_change,
        top_score_metric=args.top_score_metric,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        binding_percentile_threshold=args.binding_percentile_threshold,
        immunogenicity_percentile_threshold=args.immunogenicity_percentile_threshold,
        presentation_percentile_threshold=args.presentation_percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
    ).execute()

if __name__ == "__main__":
    main()
