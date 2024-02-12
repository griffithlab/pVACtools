import tempfile
import shutil

from pvactools.lib.identify_problematic_amino_acids import IdentifyProblematicAminoAcids
from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes, PvacfuseAggregateAllEpitopes, PvacbindAggregateAllEpitopes
from pvactools.lib.binding_filter import BindingFilter
from pvactools.lib.filter import Filter, FilterCriterion
from pvactools.lib.top_score_filter import PvacseqTopScoreFilter, PvacfuseTopScoreFilter, PvacbindTopScoreFilter
from pvactools.lib.calculate_manufacturability import CalculateManufacturability
from pvactools.lib.calculate_reference_proteome_similarity import CalculateReferenceProteomeSimilarity
from pvactools.lib.net_chop import NetChop
from pvactools.lib.netmhc_stab import NetMHCStab

class PostProcessor:
    def __init__(self, **kwargs):
        for (k,v) in kwargs.items():
           setattr(self, k, v)
        self.aggregate_report = self.input_file.replace('.tsv', '.aggregated.tsv')
        self.identify_problematic_amino_acids_fh = tempfile.NamedTemporaryFile()
        self.binding_filter_fh = tempfile.NamedTemporaryFile()
        self.coverage_filter_fh = tempfile.NamedTemporaryFile()
        self.transcript_support_level_filter_fh = tempfile.NamedTemporaryFile()
        self.top_score_filter_fh = tempfile.NamedTemporaryFile()
        self.net_chop_fh = tempfile.NamedTemporaryFile()
        self.netmhc_stab_fh = tempfile.NamedTemporaryFile()
        self.manufacturability_fh = tempfile.NamedTemporaryFile()
        self.reference_similarity_fh = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.file_type = kwargs.pop('file_type', None)
        self.fasta = kwargs.pop('fasta', None)
        self.net_chop_fasta = kwargs.pop('net_chop_fasta', None)
        if not hasattr(self, 'flurry_state'):
            self.flurry_state = self.get_flurry_state()
        self.el_only = all([self.is_el(a) for a in self.prediction_algorithms])

    def get_flurry_state(self):
        if 'MHCflurry' in self.prediction_algorithms and 'MHCflurryEL' in self.prediction_algorithms:
            self.prediction_algorithms.remove('MHCflurryEL')
            return 'both'
        elif 'MHCflurry' in self.prediction_algorithms:
            return 'BA_only'
        elif 'MHCflurryEL' in self.prediction_algorithms:
            pred_idx = self.prediction_algorithms.index('MHCflurryEL')
            self.prediction_algorithms[pred_idx] = 'MHCflurry'
            return 'EL_only'
        else:
            return None

    def is_el(self, algorithm):
        if algorithm == 'MHCflurry' and self.flurry_state == 'EL_only':
            return True
        if algorithm in ['NetMHCIIpanEL', 'NetMHCpanEL', 'BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
            return True
        return False

    def execute(self):
        self.identify_problematic_amino_acids()
        self.aggregate_all_epitopes()
        self.calculate_manufacturability()
        self.execute_binding_filter()
        self.execute_coverage_filter()
        self.execute_transcript_support_level_filter()
        self.execute_top_score_filter()
        self.call_net_chop()
        self.call_netmhc_stab()
        self.calculate_reference_proteome_similarity()
        if not self.el_only:
            shutil.copy(self.reference_similarity_fh.name, self.aggregate_report)
        shutil.copy(self.netmhc_stab_fh.name, self.filtered_report_file)
        self.close_filehandles()
        print("\nDone: Pipeline finished successfully. File {} contains list of filtered putative neoantigens.\n".format(self.filtered_report_file))

    def identify_problematic_amino_acids(self):
        if self.problematic_amino_acids:
            print("Identifying peptides with problematic amino acids")
            IdentifyProblematicAminoAcids(self.input_file, self.identify_problematic_amino_acids_fh.name, self.problematic_amino_acids, file_type=self.file_type).execute()
            shutil.copy(self.identify_problematic_amino_acids_fh.name, self.input_file)
            print("Completed")

    def aggregate_all_epitopes(self):
        if self.el_only:
            return
        print("Creating aggregated report")
        if self.file_type == 'pVACseq':
            PvacseqAggregateAllEpitopes(
                self.input_file,
                self.aggregate_report,
                tumor_purity=self.tumor_purity,
                binding_threshold=self.binding_threshold,
                percentile_threshold=self.percentile_threshold,
                allele_specific_binding_thresholds=self.allele_specific_binding_thresholds,
                trna_vaf=self.trna_vaf,
                trna_cov=self.trna_cov,
                expn_val=self.expn_val,
                maximum_transcript_support_level=self.maximum_transcript_support_level,
                top_score_metric=self.top_score_metric,
                allele_specific_anchors=self.allele_specific_anchors,
                anchor_contribution_threshold=self.anchor_contribution_threshold,
                aggregate_inclusion_binding_threshold=self.aggregate_inclusion_binding_threshold,
            ).execute()
        elif self.file_type == 'pVACfuse':
            PvacfuseAggregateAllEpitopes(
                self.input_file,
                self.aggregate_report,
                binding_threshold=self.binding_threshold,
                allele_specific_binding_thresholds=self.allele_specific_binding_thresholds,
                percentile_threshold=self.percentile_threshold,
                top_score_metric=self.top_score_metric,
                read_support=self.read_support,
                expn_val=self.expn_val,
                aggregate_inclusion_binding_threshold=self.aggregate_inclusion_binding_threshold,
            ).execute()
        elif self.file_type == 'pVACbind':
            PvacbindAggregateAllEpitopes(
                self.input_file,
                self.aggregate_report,
                binding_threshold=self.binding_threshold,
                allele_specific_binding_thresholds=self.allele_specific_binding_thresholds,
                percentile_threshold=self.percentile_threshold,
                top_score_metric=self.top_score_metric,
                aggregate_inclusion_binding_threshold=self.aggregate_inclusion_binding_threshold,
            ).execute()
        print("Completed")

    def calculate_manufacturability(self):
        if self.run_manufacturability_metrics:
            print("Calculating Manufacturability Metrics")
            CalculateManufacturability(self.input_file, self.manufacturability_fh.name, self.file_type).execute()
            shutil.copy(self.manufacturability_fh.name, self.input_file)
            print("Completed")

    def execute_binding_filter(self):
        if self.el_only:
            shutil.copy(self.input_file, self.binding_filter_fh.name)
            return
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
            if self.file_type == 'pVACseq':
                filter_criteria.append(FilterCriterion("Normal Depth", '>=', self.normal_cov, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Normal VAF", '<=', self.normal_vaf, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Tumor DNA Depth", '>=', self.tdna_cov, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Tumor DNA VAF", '>=', self.tdna_vaf, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Tumor RNA Depth", '>=', self.trna_cov, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Tumor RNA VAF", '>=', self.trna_vaf, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Gene Expression", '>=', self.expn_val, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Transcript Expression", '>=', self.expn_val, exclude_nas=self.exclude_NAs))
            elif self.file_type == 'pVACfuse':
                filter_criteria.append(FilterCriterion("Read Support", '>=', self.read_support, exclude_nas=self.exclude_NAs))
                filter_criteria.append(FilterCriterion("Expression", '>=', self.expn_val, exclude_nas=self.exclude_NAs))
            Filter(self.binding_filter_fh.name, self.coverage_filter_fh.name, filter_criteria).execute()
            print("Completed")
        else:
            shutil.copy(self.binding_filter_fh.name, self.coverage_filter_fh.name)

    def execute_transcript_support_level_filter(self):
        if self.run_transcript_support_level_filter:
            print("Running Transcript Support Level Filter")
            filter_criteria = [FilterCriterion('Transcript Support Level', '<=', self.maximum_transcript_support_level, exclude_nas=True, skip_value='Not Supported')]
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
        if self.el_only:
            shutil.copy(self.transcript_support_level_filter_fh.name, self.top_score_filter_fh.name)
            return
        print("Running Top Score Filter")
        if self.file_type == 'pVACseq':
            PvacseqTopScoreFilter(
                self.transcript_support_level_filter_fh.name,
                self.top_score_filter_fh.name,
                top_score_metric=self.top_score_metric,
                binding_threshold=self.binding_threshold,
                allele_specific_binding_thresholds=self.allele_specific_binding_thresholds,
                maximum_transcript_support_level=self.maximum_transcript_support_level,
                allele_specific_anchors=self.allele_specific_anchors,
                anchor_contribution_threshold=self.anchor_contribution_threshold,
            ).execute()
        elif self.file_type == 'pVACfuse':
            PvacfuseTopScoreFilter(
                self.transcript_support_level_filter_fh.name,
                self.top_score_filter_fh.name,
                top_score_metric = self.top_score_metric,
            ).execute()
        elif self.file_type == 'pVACbind':
            PvacbindTopScoreFilter(
                self.transcript_support_level_filter_fh.name,
                self.top_score_filter_fh.name,
                top_score_metric = self.top_score_metric,
            ).execute()
        print("Completed")

    def call_net_chop(self):
        if self.run_net_chop:
            print("Submitting remaining epitopes to NetChop")
            NetChop(self.top_score_filter_fh.name, self.net_chop_fasta, self.net_chop_fh.name, self.net_chop_method, str(self.net_chop_threshold), self.file_type).execute()
            print("Completed")
        else:
            shutil.copy(self.top_score_filter_fh.name, self.net_chop_fh.name)

    def call_netmhc_stab(self):
        if self.run_netmhc_stab:
            print("Running NetMHCStabPan")
            NetMHCStab(self.net_chop_fh.name, self.netmhc_stab_fh.name, self.file_type, self.top_score_metric).execute()
            print("Completed")
        else:
            shutil.copy(self.net_chop_fh.name, self.netmhc_stab_fh.name)

    def calculate_reference_proteome_similarity(self):
        if self.el_only:
            return
        if self.run_reference_proteome_similarity:
            print("Calculating Reference Proteome Similarity")
            if self.file_type == 'pVACseq':
                aggregate_metrics_file = self.aggregate_report.replace('.tsv', '.metrics.json')
                CalculateReferenceProteomeSimilarity(
                    self.aggregate_report,
                    self.fasta,
                    self.reference_similarity_fh.name,
                    species=self.species,
                    file_type=self.file_type,
                    n_threads=self.n_threads,
                    blastp_path=self.blastp_path,
                    blastp_db=self.blastp_db,
                    peptide_fasta=self.peptide_fasta,
                    aggregate_metrics_file=aggregate_metrics_file,
                ).execute()
                aggregate_metrics_output_file = self.reference_similarity_fh.name.replace('.tsv', '.metrics.json')
                shutil.move(aggregate_metrics_output_file, aggregate_metrics_file)
            else:
                CalculateReferenceProteomeSimilarity(
                    self.aggregate_report,
                    self.fasta,
                    self.reference_similarity_fh.name,
                    species=self.species,
                    file_type=self.file_type,
                    n_threads=self.n_threads,
                    blastp_path=self.blastp_path,
                    blastp_db=self.blastp_db,
                    peptide_fasta=self.peptide_fasta,
                ).execute()
            shutil.move("{}.reference_matches".format(self.reference_similarity_fh.name), "{}.reference_matches".format(self.aggregate_report))
            print("Completed")
        else:
            shutil.copy(self.aggregate_report, self.reference_similarity_fh.name)

    def close_filehandles(self):
        self.binding_filter_fh.close()
        self.coverage_filter_fh.close()
        self.transcript_support_level_filter_fh.close()
        self.top_score_filter_fh.close()
        self.net_chop_fh.close()
        self.netmhc_stab_fh.close()
        self.manufacturability_fh.close()
        self.reference_similarity_fh.close()
