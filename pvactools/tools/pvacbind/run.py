import sys
import argparse
import os
import shutil
import yaml
import platform

from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import PvacbindPipeline
from pvactools.lib.run_argument_parser import PvacbindRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.tools.pvacbind.run_planning import (
    CombinedReportHandoff,
    PvacbindOutputLayout,
    PvacbindPipelineHandoffBuilder,
    PvacbindRunPlan,
)

def define_parser():
    return PvacbindRunArgumentParser().parser

def create_combined_reports(base_output_dir, args):
    output_layout = PvacbindOutputLayout(base_output_dir, args.sample_name)
    combined_handoff = CombinedReportHandoff(output_layout)
    os.makedirs(output_layout.combined_output_dir, exist_ok=True)

    if not os.path.exists(output_layout.class_i_all_epitopes_file):
        print("File {} doesn't exist. Aborting.".format(output_layout.class_i_all_epitopes_file))
        return
    if not os.path.exists(output_layout.class_ii_all_epitopes_file):
        print("File {} doesn't exist. Aborting.".format(output_layout.class_ii_all_epitopes_file))
        return

    combine_reports(combined_handoff.input_files, output_layout.combined_all_epitopes_file)
    post_processing_params = combined_handoff.post_processing_params(args)

    PostProcessor(**post_processing_params).execute()

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    if args.n_threads > 1 and platform.system() == "Darwin":
        raise Exception("Multithreading is not supported on MacOS")

    base_output_dir = os.path.abspath(args.output_dir)

    if (args.netmhciipan_version == '4.0' and args.iedb_install_directory is not None):
        raise Exception("Standalone IEDB does not support version 4.0")
    NetMHCIIVersion.netmhciipan_version = args.netmhciipan_version

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    alleles = combine_class_ii_alleles(args.allele)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(alleles)
    output_layout = PvacbindOutputLayout(base_output_dir, args.sample_name)
    run_plan = PvacbindRunPlan(
        class_i_prediction_algorithms,
        class_ii_prediction_algorithms,
        class_i_alleles,
        class_ii_alleles,
        species,
    )
    pipeline_handoffs = PvacbindPipelineHandoffBuilder(args, run_plan, output_layout)

    if run_plan.should_run_class_i():
        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None
        
        if args.use_normalized_percentiles and species != 'human':
            print("WARNING: Normalized percentiles are only available for human alleles. Option will be ignored.")
            args.use_normalized_percentiles = False

        print("Executing MHC Class I predictions")

        output_dir = output_layout.class_i_output_dir
        os.makedirs(output_dir, exist_ok=True)

        class_i_arguments = pipeline_handoffs.class_i_arguments(iedb_mhc_i_executable)
        pipeline = PvacbindPipeline(**class_i_arguments)
        pipeline.execute()
    else:
        print(run_plan.class_i_skip_message())

    if run_plan.should_run_class_ii():
        if args.iedb_install_directory:
            iedb_mhc_ii_executable = os.path.join(args.iedb_install_directory, 'mhc_ii', 'mhc_II_binding.py')
            if not os.path.exists(iedb_mhc_ii_executable):
                sys.exit("IEDB MHC II executable path doesn't exist %s" % iedb_mhc_ii_executable)
        else:
            iedb_mhc_ii_executable = None

        print("Executing MHC Class II predictions")

        output_dir = output_layout.class_ii_output_dir
        os.makedirs(output_dir, exist_ok=True)

        class_ii_arguments = pipeline_handoffs.class_ii_arguments(iedb_mhc_ii_executable)
        pipeline = PvacbindPipeline(**class_ii_arguments)
        pipeline.execute()
    else:
        print(run_plan.class_ii_skip_message())

    if run_plan.should_create_combined_reports():
        print("Creating combined reports")
        create_combined_reports(base_output_dir, args)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
