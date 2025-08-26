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
    ):
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.anchor_calculator = anchor_calculator
        self.mt_top_score_metric = mt_top_score_metric
        self.top_score_mode = top_score_mode

    def get(self, df):
        #get all entries with Biotype 'protein_coding'
        biotype_df = df[df['Biotype'] == 'protein_coding']
        #if there are none, reset to previous dataframe
        if biotype_df.shape[0] == 0:
            biotype_df = df

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

        #subset prob_pos dataframe to only include entries that pass the anchor position check
        prob_pos_df['anchor_residue_pass'] = prob_pos_df.apply(lambda x: self.anchor_calculator.is_anchor_residue_pass(x), axis=1)
        anchor_residue_pass_df = prob_pos_df[prob_pos_df['anchor_residue_pass']]
        if anchor_residue_pass_df.shape[0] == 0:
            anchor_residue_pass_df = prob_pos_df

        #determine the entry with the lowest IC50 Score, transcript prioritization status, and longest Transcript
        anchor_residue_pass_df.sort_values(by=[
            "{} MT {}".format(self.mt_top_score_metric, self.top_score_mode),
            "Transcript Length",
        ], inplace=True, ascending=[True, False])
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
        if 'Expression' in df:
            df['Expression Sort'] = df['Expression']
            df['Expression Sort'].replace({'NA': 0})
        prob_pos_df.sort_values(by=["{} {}".format(self.top_score_metric, self.top_score_mode), 'Expression Sort'], inplace=True, ascending=[True, False])
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
    ):
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.mt_top_score_metric = mt_top_score_metric
        self.top_score_mode = top_score_mode

    def get(self, df):
        #get all entries with Biotype 'protein_coding'
        biotype_df = df[df['Biotype'] == 'protein_coding']
        #if there are none, reset to previous dataframe
        if biotype_df.shape[0] == 0:
            biotype_df = df

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

        #determine the entry with the lowest IC50 Score, transcript prioritization status, and longest Transcript
        prob_pos_df.sort_values(by=[
            "{} {}".format(self.mt_top_score_metric, self.top_score_mode),
            "WT Protein Length",
        ], inplace=True, ascending=[True, False])
        return prob_pos_df.iloc[0]
