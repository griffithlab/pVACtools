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
    """Estimate the maximum size a column in the data table will need"""
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
    return 10 + (8*size)

validation_cache = {} #saves time
def validate_column(iterable, column):
    """Returns True iff at least one value in the column is not None"""
    if column not in validation_cache:
        validation_cache[column] = functools.reduce(
            lambda x,y: x or y is not None,
            iterable,
        )
    return validation_cache[column]

def create_filters(filters):
    """Takes an iterable list of unary filtering operators\
    and returns a combined filter function"""
    def apply_filters(entry):
        result = functools.reduce(
            lambda x,y: x and y(entry),
            filters,
            True
        )
        if not result:
            apply_filters.count += 1
        return result
    apply_filters.count = 0
    return apply_filters

def make_filter(colname, getter):
    """Convenience method for generating min/max filters for a column"""
    return [
        lambda x:x[colname] is None or x[colname] >= getter()[0],
        lambda x:x[colname] is None or x[colname] <= getter()[1]
    ]

def stepUp(val, step):
    """Convenience method for picking the maximum range of a filter"""
    return (int(val/step)+1)/(1/step) # avoids more floating point errors than multiplying by step

def range_column_filter(colname, stepsize, title=None):
    """Create a column filter, if the column exists and contains at least one\
    not-None value.  Creates the slider for the filter and returns a getter\
    for the slider's current range"""
    if colname in entries[0]:
        column_data = [entry[colname] for entry in entries]
        if validate_column(column_data, colname):
            top = stepUp(max((val for val in column_data if val is not None)), stepsize)
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

#Parse session arguments
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
#Fetch table data from postgres
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
from bokeh.models import TableColumn, TapTool, BoxSelectTool, ResizeTool
from bokeh.models.ranges import Range1d as Range
from bokeh.models.widgets import Select, DataTable, Toggle, RadioButtonGroup
from bokeh.plotting import figure
#Set up the x/y axis selectors and the toggle to hide null entries
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
hide_null = Toggle(
    active=True,
    label="Hide 0 results with null X or Y axis values"
)
widgets.append(hide_null)

presets = RadioButtonGroup(
    labels=["MT vs WT Epitope Affinity", "Tumor Clonality and Expression"], active=0)

def available(x):
    for entry in entries:
        if x in entry and entry[x] is not None:
            return True
    return False

if not available('corresponding_wt_score') or not available('best_mt_score'):
    presets.labels.remove('MT vs WT Epitope Affinity')
if not available('tumor_dna_vaf') or not available('tumor_rna_vaf'):
    presets.labels.remove('Tumor Clonality and Expression')
widgets.append(presets)

#Set up the data dictionary (a transposed version of entries)
data_dict = {
    key:[] for key in cols
}
data_dict.update({
    '_x':[],
    '_y':[],
})
source = ColumnDataSource(data=data_dict) #wrap a datasource around the dictionary

p = figure(
    title = sample,
    # sizing_mode='stretch_both',
    plot_height=800, plot_width=900,
)
#every keyword argument can accept a constant value, or a column name
#if given a column name, it will use the values of that column in the data source
#for each point
p.circle(x="_x", y="_y", source=source, size=7, color="blue", line_color=None, fill_alpha=1)
p.add_tools(TapTool())
p.add_tools(BoxSelectTool())
p.add_tools(ResizeTool())
hover = HoverTool()
hover.tooltips = [
    ('Row', '@rowid'),
    ('X', '@_x'),
    ('Y', '@_y')
]
p.add_tools(hover)
#add the data table
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
        for colname in sorted(cols, key=lambda x:'' if x=='rowid' else x)
        if validate_column((entry[colname] for entry in entries), colname)
    ],
    fit_columns = False,
    selectable = True,
    source = source,
    row_headers = False,
    # sizing_mode = 'scale_width',
    width = 1200
)
#Update function
def update():
    """Updates the data in the datasource for the graph and table based on\
    the currently applied filters"""
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
    filtered = [
        entry for entry in entries
        if (
            not (entry[x] is None and hide_null.active) and
            not (entry[y] is None and hide_null.active) and
            data_filter(entry)
        )
    ]
    data_dict = {
        key:[
            entry[key] for entry in filtered
        ]
        for key in cols
    }
    data_dict.update({
        '_x':[entry for entry in data_dict[x]],
        '_y':[entry for entry in data_dict[y]]
    })
    hide_null.label="%s %d results with null X or Y axis values"%(
        'Show' if hide_null.active else 'Hide',
        len([entry for entry in entries if(entry[x] is None or entry[y] is None)])
    )
    source.data = data_dict
    #end update function

#create range filters for various columns
getters.append((
    'best_mt_score',
    range_column_filter('best_mt_score', 50, 'Binding Threshold (best)')
))
getters.append((
    'median_mt_score',
    range_column_filter('median_mt_score', 50, 'Binding Threshold (median)')
))
getters.append((
    'corresponding_wt_score',
    range_column_filter('corresponding_wt_score', 50, 'Binding Threshold (WT)')
))
getters.append((
    'corresponding_fold_change',
    range_column_filter('corresponding_fold_change', .1, 'Fold Change')
))
getters.append((
    'normal_depth',
    range_column_filter('normal_depth', 5)
))
getters.append((
    'normal_vaf',
    range_column_filter('normal_vaf', .5)
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
getters.append((
    'tumor_rna_vaf',
    range_column_filter('tumor_rna_vaf', .5)
))
getters.append((
    'best_cleavage_score',
    range_column_filter('best_cleavage_score', .05)
))
getters.append((
    'half_life',
    range_column_filter('half_life', .05)
))
getters.append((
    'predicted_stability',
    range_column_filter('predicted_stability', .05)
))
getters.append((
    'gene_expression',
    range_column_filter('gene_expression', 1)
))
getters.append((
    'transcript_expression',
    range_column_filter('transcript_expression', 1)
))

def change_preset():
    if presets.labels[presets.active] == "MT vs WT Epitope Affinity":
        for item in widgets:
            try:
                if item.title == 'X-Axis Value':
                    item.value = 'corresponding_wt_score'
                elif item.title == 'Y-Axis Value':
                    item.value = 'best_mt_score'
            except AttributeError:
                continue

    if presets.labels[presets.active] == "Tumor Clonality and Expression":
        for item in widgets:
            try:
                if item.title == 'X-Axis Value':
                    item.value = 'tumor_dna_vaf'
                elif item.title == 'Y-Axis Value':
                    item.value = 'tumor_rna_vaf'
            except AttributeError:
                continue

#Add callbacks to the 3 widgets manually created back at the start
x_field.on_change('value', lambda a,r,g: update())
y_field.on_change('value', lambda a,r,g: update())
presets.on_change('active', lambda a,r,g: change_preset())
hide_null.on_click(lambda arg: update())

#Add all models and widgets to the document
box = widgetbox(*widgets, sizing_mode='stretch_both')
fig = column(
    row(box, p),
    table,
    sizing_mode='scale_width'
)
update() #initial update
curdoc().add_root(fig)
curdoc().title = sample
