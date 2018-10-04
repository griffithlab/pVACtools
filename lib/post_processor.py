import tempfile
import shutil
from lib.binding_filter import *
from lib.filter import *
from lib.top_score_filter import *
from lib.condense_final_report import *
from lib.rank_epitopes import *
from lib.post_processor import *
import lib.net_chop
import lib.netmhc_stab

class PostProcessor:
    def __init__(self, **kwargs):
        for (k,v) in kwargs.items():
           setattr(self, k, v)
        self.binding_filter_fh = tempfile.NamedTemporaryFile()
        self.coverage_filter_fh = tempfile.NamedTemporaryFile()
        self.transcript_support_level_filter_fh = tempfile.NamedTemporaryFile()
        self.top_score_filter_fh = tempfile.NamedTemporaryFile()
        self.net_chop_fh = tempfile.NamedTemporaryFile()
        self.netmhc_stab_fh = tempfile.NamedTemporaryFile()
        self.condensed_report_fh = tempfile.NamedTemporaryFile()
        self.ranked_epitopes_fh = tempfile.NamedTemporaryFile()

    def execute(self):
        self.execute_binding_filter()
        self.execute_coverage_filter()
        self.execute_transcript_support_level_filter()
        self.execute_top_score_filter()
        self.call_net_chop()
        self.call_netmhc_stab()
        self.condense_report()
        self.rank_epitopes()
        shutil.copy(self.netmhc_stab_fh.name, self.filtered_report_file)
        shutil.copy(self.ranked_epitopes_fh.name, self.condensed_report_file)
        self.close_filehandles()

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
        ).execute()
        print("Completed")

    def execute_coverage_filter(self):
        if self.run_coverage_filter:
            print("Running Coverage Filters")
            filter_criteria = []
            filter_criteria.append({'column': "Normal_Depth", 'operator': '>=', 'threshold': self.normal_cov})
            filter_criteria.append({'column': "Normal_VAF", 'operator': '<=', 'threshold': self.normal_vaf})
            filter_criteria.append({'column': "Tumor_DNA_Depth", 'operator': '>=', 'threshold': self.tdna_cov})
            filter_criteria.append({'column': "Tumor_DNA_VAF", 'operator': '>=', 'threshold': self.tdna_vaf})
            filter_criteria.append({'column': "Tumor_RNA_Depth", 'operator': '>=', 'threshold': self.trna_cov})
            filter_criteria.append({'column': "Tumor_RNA_VAF", 'operator': '>=', 'threshold': self.trna_vaf})
            filter_criteria.append({'column': "Gene_Expression", 'operator': '>=', 'threshold': self.expn_val})
            filter_criteria.append({'column': "Transcript_Expression", 'operator': '>=', 'threshold': self.expn_val})
            Filter(self.binding_filter_fh.name, self.coverage_filter_fh.name, filter_criteria, self.exclude_NAs).execute()
            print("Completed")
        else:
            shutil.copy(self.binding_filter_fh.name, self.coverage_filter_fh.name)

    def execute_transcript_support_level_filter(self):
        if self.run_transcript_support_level_filter:
            print("Running Transcript Support Level Filter")
            filter_criteria = [{'column': 'Transcript Support Level', 'operator': '<=', 'threshold': self.maximum_transcript_support_level}]
            Filter(
                self.coverage_filter_fh.name,
                self.transcript_support_level_filter_fh.name,
                filter_criteria,
                self.exclude_NAs,
                ['Transcript Support Level'],
            ).execute()
            print("Complete")
        else:
            shutil.copy(self.coverage_filter_fh.name, self.transcript_support_level_filter_fh.name)

    def execute_top_score_filter(self):
        print("Running Top Score Filter")
        TopScoreFilter(self.transcript_support_level_filter_fh.name, self.top_score_filter_fh.name, self.top_score_metric).execute()
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

    def condense_report(self):
        print("Creating Condensed Report")
        CondenseFinalReport(self.netmhc_stab_fh.name, self.condensed_report_fh.name, self.top_score_metric).execute()
        print("Completed")

    def rank_epitopes(self):
        print("Ranking neoepitopes")
        RankEpitopes(self.condensed_report_fh.name, self.ranked_epitopes_fh.name).execute()
        print("Completed")

    def close_filehandles(self):
        self.binding_filter_fh.close()
        self.coverage_filter_fh.close()
        self.transcript_support_level_filter_fh.close()
        self.top_score_filter_fh.close()
        self.net_chop_fh.close()
        self.netmhc_stab_fh.close()
        self.condensed_report_fh.close()
        self.ranked_epitopes_fh.close()
