import argparse
import io
import os
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest import mock

import pvactools.tools.pvacbind.run as pvacbind_run


class PvacbindRunOrchestrationTests(unittest.TestCase):
    def test_main_builds_class_specific_pipeline_handoffs_and_combined_report_request(self):
        input_file = "/input.fa"
        pipeline_calls = []
        executed_pipelines = []
        postprocessor_calls = []

        class FakePipeline:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                pipeline_calls.append(kwargs)

            def execute(self):
                executed_pipelines.append(self.kwargs)
                all_epitopes_path = os.path.join(
                    self.kwargs["output_dir"],
                    "{}.{}.all_epitopes.tsv".format(
                        self.kwargs["sample_name"],
                        self.kwargs["filename_addition"],
                    ),
                )
                with open(all_epitopes_path, "w"):
                    pass

        class FakePostProcessor:
            def __init__(self, **kwargs):
                postprocessor_calls.append(kwargs)

            def execute(self):
                pass

        with tempfile.TemporaryDirectory() as output_dir, \
                mock.patch.object(pvacbind_run, "PvacbindPipeline", FakePipeline), \
                mock.patch.object(pvacbind_run, "combine_reports") as combine_reports, \
                mock.patch.object(pvacbind_run, "PostProcessor", FakePostProcessor), \
                mock.patch.object(pvacbind_run, "change_permissions_recursive") as change_permissions_recursive, \
                mock.patch.object(pvacbind_run.platform, "system", return_value="Linux"):
            pvacbind_run.main([
                input_file,
                "Sample",
                "HLA-A*02:01,DRB1*11:01",
                "NetMHC",
                "NNalign",
                output_dir,
                "-e1", "8,9",
                "-e2", "15,16",
                "--top-score-metric", "lowest",
                "--top-score-metric2", "ic50,presentation_percentile",
                "--binding-threshold", "123",
                "--binding-percentile-threshold", "1.5",
                "--presentation-percentile-threshold", "3.5",
                "--immunogenicity-percentile-threshold", "4.5",
                "--percentile-threshold-strategy", "exploratory",
                "--allele-specific-binding-thresholds",
                "--net-chop-method", "cterm",
                "--net-chop-threshold", "0.7",
                "--additional-report-columns", "sample_name",
                "--fasta-size", "42",
                "--iedb-retries", "7",
                "--keep-tmp-files",
                "--n-threads", "2",
                "--run-reference-proteome-similarity",
                "--blastp-path", "/blastp",
                "--blastp-db", "refseq_protein",
                "--problematic-amino-acids", "C:1,M",
                "--peptide-fasta", "peptides.fa",
                "--aggregate-inclusion-binding-threshold", "4321",
                "--aggregate-inclusion-count-limit", "9",
                "--netmhc-stab",
                "--use-normalized-percentiles",
                "--reference-scores-path", "/scores",
            ])

        self.assertEqual(len(pipeline_calls), 2)
        self.assertEqual(executed_pipelines, pipeline_calls)

        class_i_arguments = pipeline_calls[0]
        class_ii_arguments = pipeline_calls[1]

        expected_shared_arguments = {
            "input_file": input_file,
            "input_file_type": "fasta",
            "sample_name": "Sample",
            "top_score_metric": "lowest",
            "top_score_metric2": ["ic50", "presentation_percentile"],
            "binding_threshold": 123,
            "binding_percentile_threshold": 1.5,
            "immunogenicity_percentile_threshold": 4.5,
            "presentation_percentile_threshold": 3.5,
            "percentile_threshold_strategy": "exploratory",
            "allele_specific_binding_thresholds": True,
            "net_chop_fasta": input_file,
            "net_chop_method": "cterm",
            "net_chop_threshold": 0.7,
            "additional_report_columns": "sample_name",
            "fasta_size": 42,
            "iedb_retries": 7,
            "keep_tmp_files": True,
            "n_threads": 2,
            "species": "human",
            "run_reference_proteome_similarity": True,
            "blastp_path": "/blastp",
            "blastp_db": "refseq_protein",
            "problematic_amino_acids": ["C:1", "M"],
            "run_post_processor": True,
            "peptide_fasta": "peptides.fa",
            "aggregate_inclusion_binding_threshold": 4321,
            "aggregate_inclusion_count_limit": 9,
        }
        for key, expected_value in expected_shared_arguments.items():
            self.assertEqual(class_i_arguments[key], expected_value)
            self.assertEqual(class_ii_arguments[key], expected_value)

        self.assertEqual(class_i_arguments["alleles"], ["HLA-A*02:01"])
        self.assertIsNone(class_i_arguments["iedb_executable"])
        self.assertEqual(class_i_arguments["epitope_lengths"], [8, 9])
        self.assertEqual(class_i_arguments["prediction_algorithms"], ["NetMHC"])
        self.assertEqual(class_i_arguments["output_dir"], os.path.join(os.path.abspath(output_dir), "MHC_Class_I"))
        self.assertTrue(class_i_arguments["netmhc_stab"])
        self.assertEqual(class_i_arguments["filename_addition"], "MHC_I")
        self.assertTrue(class_i_arguments["use_normalized_percentiles"])
        self.assertEqual(class_i_arguments["reference_scores_path"], "/scores")

        self.assertEqual(class_ii_arguments["alleles"], ["DRB1*11:01"])
        self.assertIsNone(class_ii_arguments["iedb_executable"])
        self.assertEqual(class_ii_arguments["epitope_lengths"], [15, 16])
        self.assertEqual(class_ii_arguments["prediction_algorithms"], ["NNalign"])
        self.assertEqual(class_ii_arguments["output_dir"], os.path.join(os.path.abspath(output_dir), "MHC_Class_II"))
        self.assertFalse(class_ii_arguments["netmhc_stab"])
        self.assertEqual(class_ii_arguments["filename_addition"], "MHC_II")
        self.assertNotIn("use_normalized_percentiles", class_ii_arguments)
        self.assertNotIn("reference_scores_path", class_ii_arguments)

        combined_dir = os.path.join(os.path.abspath(output_dir), "combined")
        class_i_file = os.path.join(os.path.abspath(output_dir), "MHC_Class_I", "Sample.MHC_I.all_epitopes.tsv")
        class_ii_file = os.path.join(os.path.abspath(output_dir), "MHC_Class_II", "Sample.MHC_II.all_epitopes.tsv")
        combined_file = os.path.join(combined_dir, "Sample.Combined.all_epitopes.tsv")
        combine_reports.assert_called_once_with([class_i_file, class_ii_file], combined_file)
        self.assertEqual(len(postprocessor_calls), 1)
        self.assertEqual(postprocessor_calls[0]["input_file"], combined_file)
        self.assertEqual(
            postprocessor_calls[0]["filtered_report_file"],
            os.path.join(combined_dir, "Sample.Combined.filtered.tsv"),
        )
        self.assertEqual(postprocessor_calls[0]["filename_addition"], "Combined")
        change_permissions_recursive.assert_called_once_with(os.path.abspath(output_dir), 0o755, 0o644)

    def test_main_skips_absent_class_work_and_does_not_request_combined_reports(self):
        pipeline_calls = []

        class FakePipeline:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                pipeline_calls.append(kwargs)

            def execute(self):
                all_epitopes_path = os.path.join(
                    self.kwargs["output_dir"],
                    "{}.{}.all_epitopes.tsv".format(
                        self.kwargs["sample_name"],
                        self.kwargs["filename_addition"],
                    ),
                )
                with open(all_epitopes_path, "w"):
                    pass

        with tempfile.TemporaryDirectory() as output_dir, \
                mock.patch.object(pvacbind_run, "PvacbindPipeline", FakePipeline), \
                mock.patch.object(pvacbind_run, "combine_reports") as combine_reports, \
                mock.patch.object(pvacbind_run, "PostProcessor") as postprocessor, \
                mock.patch.object(pvacbind_run, "change_permissions_recursive"), \
                redirect_stdout(io.StringIO()) as stdout:
            pvacbind_run.main([
                "/input.fa",
                "Sample",
                "HLA-A*02:01",
                "NetMHC",
                output_dir,
            ])

        self.assertEqual(len(pipeline_calls), 1)
        self.assertEqual(pipeline_calls[0]["filename_addition"], "MHC_I")
        combine_reports.assert_not_called()
        postprocessor.assert_not_called()
        self.assertIn(
            "No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.",
            stdout.getvalue(),
        )

    def test_create_combined_reports_builds_postprocessor_handoff(self):
        postprocessor_calls = []

        class FakePostProcessor:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                postprocessor_calls.append(kwargs)

            def execute(self):
                pass

        with tempfile.TemporaryDirectory() as base_output_dir, \
                mock.patch.object(pvacbind_run, "combine_reports") as combine_reports, \
                mock.patch.object(pvacbind_run, "PostProcessor", FakePostProcessor):
            os.makedirs(os.path.join(base_output_dir, "MHC_Class_I"))
            os.makedirs(os.path.join(base_output_dir, "MHC_Class_II"))
            class_i_file = os.path.join(base_output_dir, "MHC_Class_I", "Sample.MHC_I.all_epitopes.tsv")
            class_ii_file = os.path.join(base_output_dir, "MHC_Class_II", "Sample.MHC_II.all_epitopes.tsv")
            with open(class_i_file, "w"):
                pass
            with open(class_ii_file, "w"):
                pass

            args = argparse.Namespace(
                sample_name="Sample",
                binding_threshold=500,
                run_reference_proteome_similarity=True,
            )
            pvacbind_run.create_combined_reports(base_output_dir, args)

        combined_dir = os.path.join(base_output_dir, "combined")
        combined_all_epitopes = os.path.join(combined_dir, "Sample.Combined.all_epitopes.tsv")
        combined_filtered = os.path.join(combined_dir, "Sample.Combined.filtered.tsv")
        combine_reports.assert_called_once_with([class_i_file, class_ii_file], combined_all_epitopes)

        self.assertEqual(len(postprocessor_calls), 1)
        postprocessor_kwargs = postprocessor_calls[0]
        self.assertEqual(postprocessor_kwargs["input_file"], combined_all_epitopes)
        self.assertEqual(postprocessor_kwargs["filtered_report_file"], combined_filtered)
        self.assertFalse(postprocessor_kwargs["run_coverage_filter"])
        self.assertIsNone(postprocessor_kwargs["minimum_fold_change"])
        self.assertEqual(postprocessor_kwargs["file_type"], "pVACbind")
        self.assertFalse(postprocessor_kwargs["run_transcript_support_level_filter"])
        self.assertFalse(postprocessor_kwargs["run_net_chop"])
        self.assertFalse(postprocessor_kwargs["run_netmhc_stab"])
        self.assertFalse(postprocessor_kwargs["run_manufacturability_metrics"])
        self.assertFalse(postprocessor_kwargs["run_reference_proteome_similarity"])
        self.assertEqual(postprocessor_kwargs["filename_addition"], "Combined")

    def test_create_combined_reports_aborts_when_class_report_is_missing(self):
        with tempfile.TemporaryDirectory() as base_output_dir, \
                mock.patch.object(pvacbind_run, "combine_reports") as combine_reports, \
                mock.patch.object(pvacbind_run, "PostProcessor") as postprocessor, \
                redirect_stdout(io.StringIO()) as stdout:
            os.makedirs(os.path.join(base_output_dir, "MHC_Class_I"))
            class_i_file = os.path.join(base_output_dir, "MHC_Class_I", "Sample.MHC_I.all_epitopes.tsv")
            with open(class_i_file, "w"):
                pass

            pvacbind_run.create_combined_reports(
                base_output_dir,
                argparse.Namespace(sample_name="Sample"),
            )

        missing_file = os.path.join(base_output_dir, "MHC_Class_II", "Sample.MHC_II.all_epitopes.tsv")
        self.assertIn("File {} doesn't exist. Aborting.".format(missing_file), stdout.getvalue())
        combine_reports.assert_not_called()
        postprocessor.assert_not_called()

    def test_iedb_retry_limit_exits_before_prediction_evidence_is_consumed(self):
        with tempfile.TemporaryDirectory() as output_dir, \
                mock.patch.object(pvacbind_run, "split_algorithms") as split_algorithms:
            with self.assertRaises(SystemExit) as cm:
                pvacbind_run.main([
                    "/input.fa",
                    "Sample",
                    "HLA-A*02:01",
                    "NetMHC",
                    output_dir,
                    "--iedb-retries", "101",
                ])

        self.assertEqual(str(cm.exception), "The number of IEDB retries must be less than or equal to 100")
        split_algorithms.assert_not_called()


if __name__ == "__main__":
    unittest.main()
