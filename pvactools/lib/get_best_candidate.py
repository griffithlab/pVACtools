import pandas as pd

from pvactools.lib.anchor_residue_pass import AnchorResiduePass
from pvactools.lib.run_utils import is_preferred_transcript

class PvacseqBestCandidate:
    def __init__(
        self,
        transcript_prioritization_strategy,
        maximum_transcript_support_level,
        anchor_calculator,
        mt_top_score_metric,
        top_score_mode,
        allow_incomplete_transcripts,
    ):
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.anchor_calculator = anchor_calculator
        self.mt_top_score_metric = mt_top_score_metric
        self.top_score_mode = top_score_mode
        self.allow_incomplete_transcripts = allow_incomplete_transcripts

    def get(self, df):
        #get all entries that don't have CDS Flags
        if self.allow_incomplete_transcripts:
            cds_df = df[df['Transcript CDS Flags'] == 'None']
            if cds_df.shape[0] == 0:
                cds_df = df
        else:
            cds_df = df

        #subset cds dataframe to only get entries with Biotype 'protein_coding'
        biotype_df = cds_df[cds_df['Biotype'] == 'protein_coding']
        #if there are none, reset to previous dataframe
        if biotype_df.shape[0] == 0:
            biotype_df = cds_df

        #subset protein_coding dataframe to only preferred transcripts
        biotype_df['transcript_pass'] = biotype_df.apply(lambda x: is_preferred_transcript(x, self.transcript_prioritization_strategy, self.maximum_transcript_support_level), axis=1)
        transcript_df = biotype_df[biotype_df['transcript_pass']]
        #if this results in an empty dataframe, reset to previous dataframe
        if transcript_df.shape[0] == 0:
            transcript_df = biotype_df

        #subset transcript dataframe to only include entries with no problematic positions
        if 'Problematic Positions' in transcript_df:
            prob_pos_df = transcript_df[transcript_df['Problematic Positions'] == "None"]
            #if this results in an empty dataframe, reset to previous dataframe
            if prob_pos_df.shape[0] == 0:
                prob_pos_df = transcript_df
        else:
            prob_pos_df = transcript_df

        #subset prob_pos dataframe to only include entries that pass the anchor position check
        prob_pos_df['anchor_residue_pass'] = prob_pos_df.apply(lambda x: self.anchor_calculator.is_anchor_residue_pass(x), axis=1)
        anchor_residue_pass_df = prob_pos_df[prob_pos_df['anchor_residue_pass']]
        if anchor_residue_pass_df.shape[0] == 0:
            anchor_residue_pass_df = prob_pos_df

        #set up sorting criteria
        anchor_residue_pass_df['mane_select_sort'] = anchor_residue_pass_df["MANE Select"].apply(lambda x: 1 if x else 2)
        anchor_residue_pass_df['canonical_sort'] = anchor_residue_pass_df["Canonical"].apply(lambda x: 1 if x else 2)
        anchor_residue_pass_df['tsl_sort'] = anchor_residue_pass_df["Transcript Support Level"].apply(lambda x: 6 if x in ['NA', 'Not Supported'] or pd.isna(x) else int(x))
        sort_columns = [
            f"{self.mt_top_score_metric} MT {self.top_score_mode}",
            "mane_select_sort",
            "canonical_sort",
            "tsl_sort",
            "Transcript Length",
            "Transcript Expression"
        ]
        sort_orders = [
            True,
            True,
            True,
            True,
            False,
            False
        ]

        #Sort the dataframe according to the criteria and pick the first (best) one
        anchor_residue_pass_df.sort_values(
            by=sort_columns,
            ascending=sort_orders,
            inplace=True
        )
        return anchor_residue_pass_df.iloc[0]

class PvacfuseBestCandidate:
    def __init__(self, top_score_metric, top_score_mode):
        self.top_score_metric = top_score_metric
        self.top_score_mode = top_score_mode

    def get(self, df):
        #subset dataframe to only include entries with no problematic positions
        if 'Problematic Positions' in df:
            prob_pos_df = df[df['Problematic Positions'] == "None"]
            #if this results in an empty dataframe, reset to previous dataframe
            if prob_pos_df.shape[0] == 0:
                prob_pos_df = df
        else:
            prob_pos_df = df
        if 'Expression' in prob_pos_df:
            prob_pos_df['Expression Sort'] = prob_pos_df['Expression']
            prob_pos_df['Expression Sort'].replace({'NA': 0})
            prob_pos_df.sort_values(by=["{} {}".format(self.top_score_metric, self.top_score_mode), 'Expression Sort'], inplace=True, ascending=[True, False])
        else:
            prob_pos_df.sort_values(by=["{} {}".format(self.top_score_metric, self.top_score_mode)], inplace=True, ascending=[True])
        return prob_pos_df.iloc[0]

class PvacbindBestCandidate:
    def __init__(self, top_score_metric, top_score_mode):
        self.top_score_metric = top_score_metric
        self.top_score_mode = top_score_mode

    def get(self, df):
        if 'Problematic Positions' in df:
            prob_pos_df = df[df['Problematic Positions'] == "None"]
            #if this results in an empty dataframe, reset to previous dataframe
            if prob_pos_df.shape[0] == 0:
                prob_pos_df = df
        else:
            prob_pos_df = df
        prob_pos_df.sort_values(by=["{} {}".format(self.top_score_metric, self.top_score_mode)], inplace=True, ascending=True)
        return prob_pos_df.iloc[0]

class PvacspliceBestCandidate:
    def __init__(
        self,
        transcript_prioritization_strategy,
        maximum_transcript_support_level,
        mt_top_score_metric,
        top_score_mode,
        allow_incomplete_transcripts,
    ):
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.mt_top_score_metric = mt_top_score_metric
        self.top_score_mode = top_score_mode
        self.allow_incomplete_transcripts=allow_incomplete_transcripts

    def get(self, df):
        #get all entries that don't have CDS Flags
        if self.allow_incomplete_transcripts:
            cds_df = df[df['Transcript CDS Flags'] == 'None']
            if cds_df.shape[0] == 0:
                cds_df = df
        else:
            cds_df = df

        #subset cds dataframe to only get entries with Biotype 'protein_coding'
        biotype_df = cds_df[cds_df['Biotype'] == 'protein_coding']
        #if there are none, reset to previous dataframe
        if biotype_df.shape[0] == 0:
            biotype_df = cds_df

        #subset protein_coding dataframe to only preferred transcripts
        biotype_df['transcript_pass'] = biotype_df.apply(lambda x: is_preferred_transcript(x, self.transcript_prioritization_strategy, self.maximum_transcript_support_level), axis=1)
        transcript_df = biotype_df[biotype_df['transcript_pass']]
        #if this results in an empty dataframe, reset to previous dataframe
        if transcript_df.shape[0] == 0:
            transcript_df = biotype_df

        #subset tsl dataframe to only include entries with no problematic positions
        if 'Problematic Positions' in transcript_df:
            prob_pos_df = transcript_df[transcript_df['Problematic Positions'] == "None"]
            #if this results in an empty dataframe, reset to previous dataframe
            if prob_pos_df.shape[0] == 0:
                prob_pos_df = transcript_df
        else:
            prob_pos_df = transcript_df

        #set up sorting criteria
        prob_pos_df['mane_select_sort'] = prob_pos_df["MANE Select"].apply(lambda x: 1 if x else 2)
        prob_pos_df['canonical_sort'] = prob_pos_df["Canonical"].apply(lambda x: 1 if x else 2)
        prob_pos_df['tsl_sort'] = prob_pos_df["Transcript Support Level"].apply(lambda x: 6 if x in ['NA', 'Not Supported'] or pd.isna(x) else int(x))
        sort_columns = [
            f"{self.mt_top_score_metric} {self.top_score_mode}",
            "mane_select_sort",
            "canonical_sort",
            "tsl_sort",
            "WT Protein Length",
            "Transcript Expression"
        ]
        sort_orders = [
            True,
            True,
            True,
            True,
            False,
            False
        ]

        #Sort the dataframe according to the criteria and pick the first (best) one
        prob_pos_df.sort_values(
            by=sort_columns,
            ascending=sort_orders,
            inplace=True
        )
        return prob_pos_df.iloc[0]
