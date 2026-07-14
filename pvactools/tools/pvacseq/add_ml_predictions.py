import os
import sys

from pvactools.lib.ml_predictor import (
    run_ml_predictions,
    define_add_ml_predictions_parser,
)


def define_parser():
    return define_add_ml_predictions_parser(tool='pvacseq')


def main(args_input=sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    # Default output_dir to same folder as Class I aggregated file
    output_dir = args.output_dir if args.output_dir is not None else os.path.dirname(os.path.abspath(args.class1_aggregated))

    run_ml_predictions(
        args.class1_aggregated,
        args.class1_all_epitopes,
        args.class2_aggregated,
        args.artifacts_path,
        output_dir,
        args.sample_name,
        args.ml_threshold_accept,
        args.ml_threshold_reject,
    )


if __name__ == "__main__":
    main()
