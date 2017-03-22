### Fetch the data from the postgres server
import postgresql as psql
from bokeh.io import curdoc
import re
import json
import decimal

def debounce(fn, value, wait, doc):
    cval = value()
    def delayed_job():
        if value() == cval:
            fn()
    doc.add_timeout_callback(delayed_job, wait)

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
from bokeh.layouts import row, widgetbox
from bokeh.charts import Scatter
from bokeh.models import ColumnDataSource, PanTool, HoverTool, Slider, RangeSlider
from bokeh.models.ranges import Range1d as Range
from bokeh.models.widgets import Select
from bokeh.plotting import figure

x_field = Select(
    title="X-Axis Value",
    options=sorted([
        (key, val) for (key, val) in cols.items()
    ], key = lambda x:x[1]),
    value = 'corresponding_wt_score'
)
y_field = Select(
    title = 'Y-Axis Value',
    options=sorted([
        (key, val) for (key, val) in cols.items()
    ], key = lambda x:x[1]),
    value = 'best_mt_score'
)
binding_filter = RangeSlider(
    title = "Binding Threshold",
    # value = 50000,
    range=(0, 50000),
    start = 0,
    end = max(entry['best_mt_score'] for entry in entries),
    step=50
)
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
p.circle(x="_x", y="_y", source=source, size=7, color="blue", line_color=None, fill_alpha=1)

def update():
    binding_threshold = binding_filter.range
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
    data_dict = {
        key:[
            entry[key] for entry in entries
            if (
                entry['best_mt_score'] >= binding_threshold[0] and
                entry['best_mt_score'] <= binding_threshold[1] and
                entry[x] is not None and
                entry[y] is not None
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


x_field.on_change('value', lambda a,o,n: update())
y_field.on_change('value', lambda a,o,n: update())
binding_filter.on_change(
    'range',
    lambda a,o,n: debounce(
        update,
        lambda :binding_filter.range,
        150,
        curdoc()
    )
)
box = widgetbox(
    x_field,
    y_field,
    binding_filter
)
figure = row(box, p)
# p.add_tools(HoverTool())
update()
curdoc().add_root(figure)
curdoc().title = sample
