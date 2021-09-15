import sys

from pvactools.lib.binding_filter import BindingFilter

def define_parser():
    return BindingFilter.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    BindingFilter(args.input_file, args.output_file, args.binding_threshold, None, args.top_score_metric, args.exclude_NAs, args.allele_specific_binding_thresholds, args.percentile_threshold, 'pVACbind').execute()

if __name__ == "__main__":
    main()
