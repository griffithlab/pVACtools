### Fetch the data from the postgres server
import postgresql as psql
from bokeh.io import curdoc
import re
import json
import decimal
import functools

def debounce(fn, value, wait, doc):
    cval = value()
    def delayed_job():
        if value() == cval:
            fn()
    doc.add_timeout_callback(delayed_job, wait)

def autosize(name, iterable):
    from itertools import chain
    vals = chain([name], iterable)
    size = 0
    for val in vals:
        if type(val)==str:
            cur_size = min(20, len(val.rstrip()))
        elif type(val) == int:
            cur_size = len('%d'%val)
        elif val is not None:
            cur_size = len('%0.5f'%val)
        else:
            cur_size = 1
        size = max(size, cur_size)
    return 12 + (6*size)

def validate_column(iterable):
    return functools.reduce(
        lambda x,y: x or y is not None,
        iterable,
    )

def create_filters(filters):
    def apply_filters(entry):
        return functools.reduce(
            lambda x,y: x and (entry is None or y(entry)),
            filters,
            True
        )
    return apply_filters

def make_filter(colname, getter):
    return [
        lambda x:x[colname] >= getter()[0],
        lambda x:x[colname] <= getter()[1]
    ]

def stepUp(val, step):
    return step * (int(val/step)+1)

def range_column_filter(colname, stepsize, title=None):
    if colname in entries[0]:
        column_data = [entry[colname] for entry in entries]
        if validate_column(column_data):
            top = stepUp(max(column_data), stepsize)
            col_filter = RangeSlider(
                title = cols[colname] if title is None else title,
                range=(0, top),
                start = 0,
                end = top,
                step=stepsize
            )
            getter = lambda :col_filter.range
            col_filter.on_change(
                'range',
                lambda a,r,g: debounce(
                    update,
                    getter,
                    150,
                    curdoc()
                )
            )
            widgets.append(col_filter)
            return getter

args = curdoc().session_context.request.arguments
try:
    parentID = int(args.get('target-process')[0])
    fileID = int(args.get('target-file')[0])
    cols = json.loads(args.get('cols')[0].decode())
    cols['rowid']='Row'
    sample = str(args.get('samplename')[0])
except BaseException as e:
    raise ValueError("Unable to parse the requried arguments") from e
tablekey = "data_%s_%s" % (
    (parentID if parentID >= 0 else 'dropbox'),
    fileID
)
db = psql.open('localhost/pvacseq')
raw_data = db.prepare("SELECT %s FROM %s"%(','.join(cols), tablekey))()
entries = [
    {
        col:float(val) if isinstance(val, decimal.Decimal) else val
        for (col, val) in zip(cols, entry)
    }
    for entry in raw_data
]
entries.sort(key=lambda x:x['rowid'])
del raw_data
### The data is stored in entries.
# Entries is a list of dicts, where each dict maps a cloumn name to a value
# cols is a dict mapping the column names to a display name
# sample is the sample name of the requested file
### From here to the bottom, the code can be changed to modify the plotted data
from bokeh.layouts import row, widgetbox, column
from bokeh.charts import Scatter
from bokeh.models import ColumnDataSource, PanTool, HoverTool, Slider, RangeSlider
from bokeh.models import TableColumn, TapTool
from bokeh.models.ranges import Range1d as Range
from bokeh.models.widgets import Select, DataTable
from bokeh.plotting import figure
widgets = []
getters = []
x_field = Select(
    title="X-Axis Value",
    options=sorted([
        (key, val) for (key, val) in cols.items()
    ], key = lambda x:x[1]),
    value = 'corresponding_wt_score'
)
widgets.append(x_field)
y_field = Select(
    title = 'Y-Axis Value',
    options=sorted([
        (key, val) for (key, val) in cols.items()
    ], key = lambda x:x[1]),
    value = 'best_mt_score'
)
widgets.append(y_field)
data_dict = {
    key:[] for key in cols
}
data_dict.update({
    '_x':[],
    '_y':[],
})
source = ColumnDataSource(data=data_dict)
p = figure(
    title = sample,
    plot_height=600, plot_width=800,
)
#every keyword argument can accept a constant value, or a column name
#if given a column name, it will use the values of that column in the data source
#for each point
p.circle(x="_x", y="_y", source=source, size=7, color="blue", line_color=None, fill_alpha=1)
p.add_tools(TapTool())

table = DataTable(
    columns = [
        TableColumn(
            title=cols[colname],
            field=colname,
            width = autosize(
                cols[colname],
                (entry[colname] for entry in entries)
            )
        )
        for colname in sorted(cols)
        if validate_column((entry[colname] for entry in entries))
    ],
    fit_columns = False,
    selectable = True,
    source = source,
    # sizing_mode = 'scale_width',
    width = 1200
)
def update():
    filters = []
    for (colname, getter) in getters:
        if getter:
            filters += make_filter(colname, getter)
    x = x_field.value
    y = y_field.value
    xlabel = cols[x]
    ylabel = cols[y]
    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel
    p.title.text = "%s vs %s"%(
        ylabel,
        xlabel
    )
    data_filter = create_filters(filters)
    data_dict = {
        key:[
            entry[key] for entry in entries
            if (
                entry[x] is not None and
                entry[y] is not None and
                data_filter(entry)
            )
        ]
        for key in cols
    }
    data_dict.update({
        '_x':[entry for entry in data_dict[x]],
        '_y':[entry for entry in data_dict[y]]
    })
    source.data = data_dict
    if not len(p.select(type=HoverTool)):
        hover = HoverTool()
        hover.tooltips = [
            ('Row', '@rowid'),
            ('X', '@_x'),
            ('Y', '@_y')
        ]
        p.add_tools(hover)

getters.append((
    'best_mt_score',
    range_column_filter('best_mt_score', 50, 'Binding Threshold')
))
getters.append((
    'tumor_dna_depth',
    range_column_filter('tumor_dna_depth', 5)
))
getters.append((
    'tumor_dna_vaf',
    range_column_filter('tumor_dna_vaf', .5)
))
getters.append((
    'tumor_rna_depth',
    range_column_filter('tumor_rna_depth', 5)
))

x_field.on_change('value', lambda a,r,g: update())
y_field.on_change('value', lambda a,r,g: update())
box = widgetbox(*widgets, sizing_mode='stretch_both')
figure = column(
    row(box, p),
    table,
    sizing_mode='scale_width'
)
update() #initial update
curdoc().add_root(figure)
curdoc().title = sample
