var data = source.data;
var samplename = source.tags[0];
var selected_rows = source.selected['1d'].indices;
// filter bokeh vars from columns (it adds _x and _y positions of points to the data)
var cols = source.column_names.filter(function(col) { return col === '_x' || col === '_y' ? false : true; });
var l = data[cols[0]].length; // grab the first col of data and count elements to get number of rows
var csv = cols.join(',') + '\n';

if(selected_rows.length === 0) {
    alert('No rows selected.');
    return;
}

for (var i = 0; i < l; i++) {
    var row = [];
    cols.map(function(col) {
        if(selected_rows.indexOf(i) > -1) {
            row.push(data[col][i]);
        }
    });

    if(row.length > 0) {
        csv = csv.concat(row.join(',')) + '\n';
    }
}

var filename = samplename + '_selected_' + new Date().toISOString() + '.csv';
var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
} else {
    var link = document.createElement("a");
    link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    link.target = "_blank";
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'));
}
