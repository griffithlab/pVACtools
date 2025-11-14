import pandas as pd
from pvactools.lib.run_utils import *

def pvacseq_sort(rows, top_score_metric, top_score_metric2, file_type='full'):
    if isinstance(rows, list):
        rows = pd.DataFrame.from_dict(rows)
    if len(rows) == 0:
        return rows

    if file_type == 'aggregated':
        tier_sorter = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        rows["rank_tier"] = rows['Tier'].map(sorter_index)
        sort_columns = ["rank_tier", "rank", "Gene", "AA Change"]
        expression_column = 'Allele Expr'
    elif file_type == 'full':
        sort_columns = ["rank", "Gene Name", "Mutation"]
        expression_column = 'Gene Expression'

    rows["rank_expr"] = pd.to_numeric(rows[expression_column], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
    rows["rank"] = rows["rank_expr"]

    if file_type == 'aggregated':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metric2_to_aggregate_column(metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]
    elif file_type == 'full':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metrics_to_column('pvacseq', top_score_metric, metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]

    rows.sort_values(by=sort_columns, inplace=True, ascending=True)

    rows.drop(labels='rank', axis=1, inplace=True)
    if file_type == 'aggregated':
        rows.drop(labels='rank_tier', axis=1, inplace=True)
    rows.drop(labels='rank_expr', axis=1, inplace=True)
    for metric2 in top_score_metric2:
        rows.drop(labels=f"rank_{metric2}", axis=1, inplace=True)

    return rows

def pvacfuse_sort(rows, top_score_metric, top_score_metric2, file_type='full'):
    if isinstance(rows, list):
        rows = pd.DataFrame.from_dict(rows)
    if len(rows) == 0:
        return rows

    if file_type == 'aggregated':
        tier_sorter = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "LowReadSupport", "LowExpr", "ProbPos", "Poor"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        rows["rank_tier"] = rows['Tier'].map(sorter_index)
        sort_columns = ["rank_tier", "rank", "ID"]
        expression_column = 'Expr'
    elif file_type == 'full':
        sort_columns = ["rank", "Mutation"]
        expression_column = 'Expression'

    rows["rank_expr"] = pd.to_numeric(rows[expression_column], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
    rows["rank"] = rows["rank_expr"]

    if file_type == 'aggregated':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metric2_to_aggregate_column(metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]
    elif file_type == 'full':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metrics_to_column('pvacfuse', top_score_metric, metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]

    rows.sort_values(by=sort_columns, inplace=True, ascending=True)

    rows.drop(labels='rank', axis=1, inplace=True)
    if file_type == 'aggregated':
        rows.drop(labels='rank_tier', axis=1, inplace=True)
    rows.drop(labels='rank_expr', axis=1, inplace=True)
    for metric2 in top_score_metric2:
        rows.drop(labels=f"rank_{metric2}", axis=1, inplace=True)

    return rows

def pvacsplice_sort(rows, top_score_metric, top_score_metric2, file_type='full'):
    if isinstance(rows, list):
        rows = pd.DataFrame.from_dict(rows)
    if len(rows) == 0:
        return rows

    if file_type == 'aggregated':
        tier_sorter = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", "LowExpr", "Subclonal", "ProbPos", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        rows["rank_tier"] = rows['Tier'].map(sorter_index)
        sort_columns = ["rank_tier", "rank", "Gene", "Transcript", "AA Change"]
        expression_column = 'Allele Expr'
    elif file_type == 'full':
        sort_columns = ["rank", "Gene Name", "Transcript", "Amino Acid Change"]
        expression_column = 'Gene Expression'

    rows["rank_expr"] = pd.to_numeric(rows[expression_column], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
    rows["rank"] = rows["rank_expr"]

    if file_type == 'aggregated':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metric2_to_aggregate_column(metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]
    elif file_type == 'full':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metrics_to_column('pvacsplice', top_score_metric, metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]

    rows.sort_values(by=sort_columns, inplace=True, ascending=True)

    rows.drop(labels='rank', axis=1, inplace=True)
    if file_type == 'aggregated':
        rows.drop(labels='rank_tier', axis=1, inplace=True)
    rows.drop(labels='rank_expr', axis=1, inplace=True)
    for metric2 in top_score_metric2:
        rows.drop(labels=f"rank_{metric2}", axis=1, inplace=True)

    return rows

def pvacbind_sort(rows, top_score_metric, top_score_metric2, file_type='full'):
    if isinstance(rows, list):
        rows = pd.DataFrame.from_dict(rows)
    if len(rows) == 0:
        return rows

    if file_type == 'aggregated':
        tier_sorter = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "ProbPos", "Poor"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        rows["rank_tier"] = rows['Tier'].map(sorter_index)
        sort_columns = ["rank_tier", "rank", "ID"]
    elif file_type == 'full':
        sort_columns = ["rank", "Mutation"]

    rows["rank"] = 0
    if file_type == 'aggregated':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metric2_to_aggregate_column(metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]
    elif file_type == 'full':
        for metric2 in top_score_metric2:
            rows[f"rank_{metric2}"] = pd.to_numeric(rows[metrics_to_column('pvacbind', top_score_metric, metric2)], errors='coerce').rank(ascending=True, method='dense', na_option='bottom')
            rows["rank"] += rows[f"rank_{metric2}"]

    rows.sort_values(by=sort_columns, inplace=True, ascending=True)

    rows.drop(labels='rank', axis=1, inplace=True)
    if file_type == 'aggregated':
        rows.drop(labels='rank_tier', axis=1, inplace=True)
    for metric2 in top_score_metric2:
        rows.drop(labels=f"rank_{metric2}", axis=1, inplace=True)

    return rows
