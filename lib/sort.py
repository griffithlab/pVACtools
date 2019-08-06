def default_sort(rows, top_score_metric):
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: (float(row['Best MT Score'])))
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Corresponding Fold Change']) if row['Corresponding Fold Change'].isdigit() else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Median MT Score'])))
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: (float(row['Median MT Score'])))
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Corresponding Fold Change']) if row['Corresponding Fold Change'].isdigit() else float('inf')), reverse=True)
        sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Best MT Score'])))
    return sorted_rows

def pvacbind_sort(rows, top_score_metric):
    if top_score_metric == 'median':
        sorted_rows = sorted(rows, key=lambda row: ( float(row['Median Score']), int(row['Sub-peptide Position']), float(row['Best Score']), row['Mutation'] ) )
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(rows, key=lambda row: ( float(row['Best Score']), int(row['Sub-peptide Position']), float(row['Median Score']), row['Mutation'] ) )
    return sorted_rows
