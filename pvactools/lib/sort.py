def default_sort(rows, top_score_metric, top_score_metric2): # Adding a comment to make this file show up as edited
    top_score_mod = "IC50 Score"
    if top_score_metric2 == "percentile":
        top_score_mod = "Percentile"
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: (float(row[f'Best MT {top_score_mod}'])))
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Corresponding Fold Change']) if str(row['Corresponding Fold Change']).isdigit() else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row[f'Median MT {top_score_mod}'])))
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: (float(row[f'Median MT {top_score_mod}'])))
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Corresponding Fold Change']) if row['Corresponding Fold Change'].isdigit() else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row[f'Best MT {top_score_mod}'])))
    return sorted_rows

def default_sort_from_pd_dict(rows, top_score_metric, top_score_metric2):
    top_score_mod = "IC50 Score"
    if top_score_metric2 == "percentile":
        top_score_mod = "Percentile"
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: row[f'Best MT {top_score_mod}'])
        sorted_rows = sorted(sorted_rows, key=lambda row: (row['Corresponding Fold Change'] if isinstance(row['Corresponding Fold Change'], float) else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: row[f'Median MT {top_score_mod}'])
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: row[f'Median MT {top_score_mod}'])
        sorted_rows = sorted(sorted_rows, key=lambda row: (row['Corresponding Fold Change'] if isinstance(row['Corresponding Fold Change'], float) else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: row[f'Best MT {top_score_mod}'])
    return sorted_rows

def pvacbind_sort(rows, top_score_metric, top_score_metric2):
    top_score_mod = "IC50 Score"
    if top_score_metric2 == "percentile":
        top_score_mod = "Percentile"
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: ( float(row[f'Median {top_score_mod}']), int(row['Sub-peptide Position']), float(row[f'Best {top_score_mod}']), row['Mutation'] ) )
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: ( float(row[f'Best {top_score_mod}']), int(row['Sub-peptide Position']), float(row[f'Median {top_score_mod}']), row['Mutation'] ) )
    return sorted_rows

def pvacsplice_sort(rows, top_score_metric, top_score_metric2):
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: ( float(row[f'Median {top_score_mod}']), float(row[f'Best {top_score_mod}']), row['Index'] ) )
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: ( float(row[f'Best {top_score_mod}']), float(row[f'Median {top_score_mod}']), row['Index'] ) )
    return sorted_rows
