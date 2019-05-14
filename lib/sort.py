def default_sort(rows, top_score_metric):
    sorted_rows = sorted(rows, key=lambda row: (int(row['Sub-peptide Position'])))
    sorted_rows = sorted(sorted_rows, key=lambda row: (float(row['Corresponding Fold Change']) if row['Corresponding Fold Change'].isdigit() else float('inf')), reverse=True)
    if top_score_metric == 'median':
        sorted_rows = sorted(
            sorted_rows,
            key=lambda row: (
                row['Gene Name'],
                row['Mutation'],
                float(row['Median MT Score']),
            )
        )
    elif top_score_metric == 'lowest':
        sorted_rows = sorted(
            sorted_rows,
            key=lambda row: (
                row['Gene Name'],
                row['Mutation'],
                float(row['Best MT Score']),
            )
        )
    return sorted_rows
