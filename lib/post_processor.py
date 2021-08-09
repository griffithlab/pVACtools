import tempfile
import shutil
from lib.aggregate_all_epitopes import *
from lib.binding_filter import *
from lib.filter import *
from lib.top_score_filter import *
from lib.calculate_manufacturability import *
from lib.calculate_reference_proteome_similarity import *
import lib.net_chop
import lib.netmhc_stab

class PostProcessor:
    def __init__(self, **kwargs):
        for (k,v) in kwargs.items():
           setattr(self, k, v)
        self.aggregate_report = self.input_file.replace('.tsv', '.aggregated.tsv')
        self.binding_filter_fh = tempfile.NamedTemporaryFile()
        self.coverage_filter_fh = tempfile.NamedTemporaryFile()
        self.transcript_support_level_filter_fh = tempfile.NamedTemporaryFile()
        self.top_score_filter_fh = tempfile.NamedTemporaryFile()
        self.net_chop_fh = tempfile.NamedTemporaryFile()
        self.netmhc_stab_fh = tempfile.NamedTemporaryFile()
        self.manufacturability_fh = tempfile.NamedTemporaryFile()
        self.reference_similarity_fh = tempfile.NamedTemporaryFile()
        self.file_type = kwargs.pop('file_type', None)
        self.fasta = kwargs.pop('fasta', None)

    def execute(self):
        self.aggregate_all_epitopes()
        self.calculate_manufacturability()
        self.execute_binding_filter()
        self.execute_coverage_filter()
        self.execute_transcript_support_level_filter()
        self.execute_top_score_filter()
        self.call_net_chop()
        self.call_netmhc_stab()
        self.calculate_reference_proteome_similarity()
        shutil.copy(self.reference_similarity_fh.name, self.filtered_report_file)
        self.close_filehandles()
        print("\nDone: Pipeline finished successfully. File {} contains list of filtered putative neoantigens.\n".format(self.filtered_report_file))

    def aggregate_all_epitopes(self):
        print("Creating aggregated report")
        AggregateAllEpitopes(self.input_file, self.aggregate_report, self.file_type).execute()
        print("Completed")

    def calculate_manufacturability(self):
        if self.run_manufacturability_metrics:
            print("Calculating Manufacturability Metrics")
            CalculateManufacturability(self.input_file, self.manufacturability_fh.name, self.file_type).execute()
            shutil.copy(self.manufacturability_fh.name, self.input_file)
            print("Completed")

    def execute_binding_filter(self):
        print("Running Binding Filters")
        BindingFilter(
            self.input_file,
            self.binding_filter_fh.name,
            self.binding_threshold,
            self.minimum_fold_change,
            self.top_score_metric,
            self.exclude_NAs,
            self.allele_specific_binding_thresholds,
            self.percentile_threshold,
            self.file_type,
        ).execute()
        print("Completed")

    def execute_coverage_filter(self):
        if self.run_coverage_filter:
            print("Running Coverage Filters")
            filter_criteria = []
            filter_criteria.append({'column': "Normal_Depth", 'operator': '>=', 'threshold': self.normal_cov, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Normal_VAF", 'operator': '<=', 'threshold': self.normal_vaf, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Tumor_DNA_Depth", 'operator': '>=', 'threshold': self.tdna_cov, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Tumor_DNA_VAF", 'operator': '>=', 'threshold': self.tdna_vaf, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Tumor_RNA_Depth", 'operator': '>=', 'threshold': self.trna_cov, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Tumor_RNA_VAF", 'operator': '>=', 'threshold': self.trna_vaf, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Gene_Expression", 'operator': '>=', 'threshold': self.expn_val, 'exclude_nas': self.exclude_NAs})
            filter_criteria.append({'column': "Transcript_Expression", 'operator': '>=', 'threshold': self.expn_val, 'exclude_nas': self.exclude_NAs})
            Filter(self.binding_filter_fh.name, self.coverage_filter_fh.name, filter_criteria).execute()
            print("Completed")
        else:
            shutil.copy(self.binding_filter_fh.name, self.coverage_filter_fh.name)

    def execute_transcript_support_level_filter(self):
        if self.run_transcript_support_level_filter:
            print("Running Transcript Support Level Filter")
            filter_criteria = [{'column': 'Transcript Support Level', 'operator': '<=', 'threshold': self.maximum_transcript_support_level, 'exclude_nas': self.exclude_NAs}]
            Filter(
                self.coverage_filter_fh.name,
                self.transcript_support_level_filter_fh.name,
                filter_criteria,
                ['Transcript Support Level'],
            ).execute()
            print("Complete")
        else:
            shutil.copy(self.coverage_filter_fh.name, self.transcript_support_level_filter_fh.name)

    def execute_top_score_filter(self):
        print("Running Top Score Filter")
        TopScoreFilter(self.transcript_support_level_filter_fh.name, self.top_score_filter_fh.name, self.top_score_metric, self.file_type).execute()
        print("Completed")

    def call_net_chop(self):
        if self.run_net_chop:
            print("Submitting remaining epitopes to NetChop")
            lib.net_chop.main([
                self.top_score_filter_fh.name,
                self.net_chop_fh.name,
                '--method',
                self.net_chop_method,
                '--threshold',
                str(self.net_chop_threshold)
            ])
            print("Completed")
        else:
            shutil.copy(self.top_score_filter_fh.name, self.net_chop_fh.name)

    def call_netmhc_stab(self):
        if self.run_netmhc_stab:
            print("Running NetMHCStabPan")
            lib.netmhc_stab.main([
                self.net_chop_fh.name,
                self.netmhc_stab_fh.name,
            ])
            print("Completed")
        else:
            shutil.copy(self.net_chop_fh.name, self.netmhc_stab_fh.name)

    def calculate_reference_proteome_similarity(self):
        if self.run_reference_proteome_similarity:
            print("Calculating Reference Proteome Similarity")
            CalculateReferenceProteomeSimilarity(
                self.netmhc_stab_fh.name,
                self.fasta,
                self.reference_similarity_fh.name,
                species=self.species,
                file_type=self.file_type
            ).execute()
            shutil.move("{}.reference_matches".format(self.reference_similarity_fh.name), "{}.reference_matches".format(self.filtered_report_file))
            print("Completed")
        else:
            shutil.copy(self.netmhc_stab_fh.name, self.reference_similarity_fh.name)

    def close_filehandles(self):
        self.binding_filter_fh.close()
        self.coverage_filter_fh.close()
        self.transcript_support_level_filter_fh.close()
        self.top_score_filter_fh.close()
        self.net_chop_fh.close()
        self.netmhc_stab_fh.close()
        self.manufacturability_fh.close()
        self.reference_similarity_fh.close()
