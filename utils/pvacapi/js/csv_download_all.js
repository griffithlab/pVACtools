var data = source.data;
var samplename = source.tags[0];
// filter bokeh vars from columns (it adds _x and _y positions of points to the data)
var cols = source.column_names.filter(function(col) { return col === '_x' || col === '_y' ? false : true; });
var l = data[cols[0]].length; // grab the first col of data and count elements to get number of rows
var csv = cols.join(',') + '\n';

for (var i = 0; i < l; i++) {
    var row = [];
    cols.map(function(col) {
        row.push(data[col][i]);
    });

    csv = csv.concat(row.join(',')) + '\n';
}

var filename = samplename + '_full_' + new Date().toISOString() + '.csv';
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
