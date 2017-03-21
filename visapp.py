import postgresql as psql
# from bokeh.io import curdoc
# req = curdoc().session_context.request
# import pdb; pdb.set_trace()
db = psql.open('localhost/pvacseq')
cols = [
    col for (col,) in
    db.prepare('SELECT column_name FROM information_schema.columns WHERE table_name = $1')('data_4_1')
]
raw_data = db.prepare("SELECT %s FROM %s"%(','.join(cols), 'data_4_1'))()
entries = [
    {col:val for (col, val) in zip(cols, entry)}
    for entry in raw_data
]

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.charts import Scatter
from bokeh.models.ranges import Range1d as Range
x = [entry['start'] for entry in entries]
y = [entry['best_mt_score'] for entry in entries]
dups = {}
for i in range(len(x)):
    key = str(x[i])+':'+str(float(y[i]))
    if key in dups:
        dups[key] += 1
    else:
        dups[key] = 1
print("Finished dups")
data = {(key,float(val)) for (key, val) in zip(x,y)}
data2 = {
    'start':[item[0] for item in data],
    'score':[item[1] for item in data],
    'size':[2*dups[str(item[0])+':'+str(float(item[1]))] for item in data]
}
print("Finished data comp")
print("Plotting %d of %d points"%(len(data), len(x)))
p = Scatter(
    data2,
    x='start',
    y='score',
    # size = 'size',
    # marker='size',
    # color='size',
    legend='top_right',
    plot_height=800, plot_width=800, title="test",
    x_range =Range(int(min(x)*.9), int(max(x)*1.1)),
    y_range = Range(int(float(min(y))*.9), int(float(max(y))*1.1)),
)
# for point in data:
#     key = str(point[0])+':'+str(float(point[1]))
#     p.scatter(
#             point[0],
#             point[1],
#             marker='circle',
#             size = 2*dups[key],
#             fill_color='blue'
#         )
# thresh = len(x)*.01
# L = len(x)
# for i in range(len(x)):
#     if(i%thresh == 0):
#         print("Status: %0.1f%%"%(i/len(x)))
#     p.scatter(
#         x[i],
#         float(y[i]),
#         marker='circle',
#         size = 5,
#         fill_color='blue'
#     )
curdoc().add_root(row(p))
