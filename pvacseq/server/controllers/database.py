import os
import re
import csv
import sys
import json
from flask import current_app
from .processes import fetch_process, is_running, gen_files_list
from .utils import initialize, savedata, column_filter

_f = re.compile(r'^\d*\.\d+$')
_i = re.compile(r'^\d+$')
queryfilters = re.compile(r'(.+)(<=?|>=?|!=|==)(.+)')

def init_column_mapping(row, schema):
    """Generate initial estimates of column data types"""
    defs = {column_filter(col):'text' for col in row}
    defs.update(schema)
    for (col, val) in row.items():
        col = column_filter(col)
        if col not in schema:
            if _f.match(val):
                try:
                    float(val)
                    print("Assigning float to",col,"based on",val)
                    defs[col] = 'decimal'
                except ValueError:
                    print("ERROR: Float mismatch:", val)
            elif _i.match(val):
                try:
                    int(val)
                    print("Assigning int to",col,"based on",val)
                    defs[col] = 'integer'
                except ValueError:
                    print("ERROR: Int mismatch:", val)
    mapping = {}
    for (col, val) in defs.items():
        if 'int' in val:
            mapping[col] = int
        elif val == 'decimal':
            mapping[col] = float
        else:
            mapping[col] = str
    return (mapping, defs)

def column_mapping(row, mapping, schema):
    """Apply filtering to the current row.
    Detect if column data types need to be changed"""
    output = {}
    changes = {}
    for (col, val) in row.items():
        col = column_filter(col)
        if col not in schema and mapping[col]==str:
            if _f.match(val):
                try:
                    float(val)
                    print("Assigning float to",col,"based on",val)
                    mapping[col]=float
                    changes[col] = float
                except ValueError:
                    print("ERROR: Float mismatch:", val)
            elif _i.match(val):
                try:
                    int(val)
                    print("Assigning int to",col,"based on",val)
                    mapping[col]=int
                    changes[col] = int
                except ValueError:
                    print("ERROR: Int mismatch:", val)
        try:
            output[col] = mapping[col](val)
        except ValueError:
            output[col] = None
    return (mapping, output, changes)


def filterfile(parentID, fileID, count, page, filters, sort, direction):
    """Gets the file ID belonging to the parent.\
    For result files, the parentID is the process ID that spawned them.\
    For dropbox files, the parentID is -1"""
    data = initialize()

    #first, generate the key
    tablekey = "data_%s_%s"%(
        (parentID if parentID >=0 else 'dropbox'),
        fileID
    )

    #check if the table exists:
    db = current_app.config['storage']['db']
    query = db.prepare("SELECT 1 FROM information_schema.tables WHERE table_name = $1")
    if not len(query(tablekey)): #table does not exist
        #Open a reader to cache the file in the database
        if parentID != -1:
            process = fetch_process(parentID, data, current_app.config['storage']['children'])
            if not process[0]:
                return (
                    {
                        "code": 400,
                        "message": "The requested process (%d) does not exist"%parentID,
                        "fields": "parentID"
                    },400
                )
            if is_running(process):
                return []
            data = gen_files_list(parentID, data)
            if fileID not in range(len(process[0]['files'])):
                return (
                    {
                        "code": 400,
                        "message": "The requested fileID (%d) does not exist for this process (%d)" %(fileID, parentID),
                        "fields": "fileID"
                    },400
                )
            raw_reader = open(process[0]['files'][fileID])
        else:
            if str(fileID) not in data['dropbox']:
                return (
                    {
                        "code": 400,
                        "message": "The requested fileID (%d) does not exist in the dropbox"%fileID,
                        "fields":"fileID"
                    }
                )
            raw_reader = open(os.path.join(
                os.path.abspath(current_app.config['files']['dropbox-dir']),
                data['dropbox'][str(fileID)]
            ))
        reader = csv.DictReader(raw_reader, delimiter='\t')

        tmp_reader = open(raw_reader.name)
        tmp = csv.DictReader(tmp_reader, delimiter='\t')
        init = next(tmp)
        tmp_reader.close()

        #Get an initial estimate of column datatypes from the first row
        (mapping, column_names) = init_column_mapping(init, current_app.config['schema'])
        tablecolumns = "\n".join( #use the estimated types to create the table
            "%s %s,"%(colname, column_names[colname])
            for colname in column_names
        )[:-1]
        CREATE_TABLE = "CREATE TABLE %s (\
            rowid SERIAL PRIMARY KEY NOT NULL,\
            %s\
        )"%(tablekey, tablecolumns)
        db.execute(CREATE_TABLE)
        #mark the table for deletion when the server shuts down
        if 'db-clean' not in current_app.config:
            current_app.config['db-clean'] = [tablekey]
        else:
            current_app.config['db-clean'].append(tablekey)
        #prepare the insertion query
        insert = db.prepare("INSERT INTO %s (%s) VALUES (%s)" %(
            tablekey,
            ','.join(column_names),
            ','.join('$%d'%i for (_,i) in zip(column_names, range(1,sys.maxsize)))
        ))
        update = "ALTER TABLE %s "%tablekey
        for row in reader:
            #process each row
            #We format the data in the row and update column data types, if necessary
            (mapping, formatted, changes) = column_mapping(row, mapping, current_app.config['schema'])
            alter_cols = []
            for (k,v) in changes.items():
                #if there were any changes to the data type, update the table
                #since we only ever update a text column to int/decimal, then
                #it's okay to nullify the data
                typ = ''
                if v == int:
                    typ = 'bigint' if k in {'start', 'stop'} else 'integer'
                elif v == float:
                    typ = 'decimal'
                alter_cols.append(
                    "ALTER COLUMN %s SET DATA TYPE %s USING null"%(
                        k,
                        typ
                    )
                )
            if len(changes):
                #Re-generate the insert statement since the data types changed
                print("Alter:",update+','.join(alter_cols))
                db.execute(update+','.join(alter_cols))
                insert = db.prepare("INSERT INTO %s (%s) VALUES (%s)" %(
                    tablekey,
                    ','.join(column_names),
                    ','.join('$%d'%i for (_,i) in zip(column_names, range(1,sys.maxsize)))
                ))
            #insert the row
            insert(*[formatted[column] for column in column_names])
        raw_reader.close()
    typequery = db.prepare("SELECT column_name, data_type FROM information_schema.columns WHERE table_name = $1")
    column_defs = typequery(tablekey)
    column_maps = {}
    for (col, typ) in column_defs:
        if 'int' in typ:
            column_maps[col] = int
        elif typ == 'numeric'or typ == 'decimal':
            column_maps[col] = float
        else:
            column_maps[col] = str
    formatted_filters = []
    for i in range(len(filters)):
        f = filters[i].strip()
        if not len(f):
            continue
        result = queryfilters.match(f)
        if not result:
            return {
                "code":400,
                "message":"Encountered an invalid filter (%s)"%f,
                "fields":"filters"
            }
        colname = column_filter(result.group(1))
        if colname not in column_maps:
            return {
                "code":400,
                "message":"Unknown column name %s"%result.group(1),
                "fields":"filters"
            }
        op = result.group(2)
        typ = column_maps[colname]
        val = None
        try:
            val = column_maps[colname](result.group(3))
        except ValueError:
            return {
                "code":400,
                "message":"Value %s cannot be formatted to match the type of column %s (%s)"%(
                    result.group(3),
                    result.group(1),
                    typ
                )
            }
        if typ == str and (op in {'==', '!='}):
            formatted_filters.append(
                json.dumps(colname) + (' not ' if '!' in op else ' ') + "LIKE '%s'"%(
                    json.dumps(val)[1:-1]
                )
            )
        else: #type is numerical
            op = op.replace('==', '=')
            formatted_filters.append(
                '%s %s %s'%(
                    json.dumps(colname),
                    op,
                    json.dumps(val)
                )
            )
    raw_query = "SELECT %s FROM %s"%(
        ','.join([k[0] for k in column_defs]),
        tablekey
    )
    if len(formatted_filters):
        raw_query += " WHERE "+" AND ".join(formatted_filters)
    if sort:
        if column_filter(sort) not in column_maps:
            return {
                'code':400,
                'message':'Invalid column name %s'%sort,
                'fields':'sort'
            }
        raw_query += " ORDER BY %s"%(column_filter(sort))
        if direction:
            raw_query += " "+direction
    if count:
        raw_query += " LIMIT %d"%count
    if page:
        raw_query += " OFFSET %d"%(page*count)
    print("Query:",raw_query)
    query = db.prepare(raw_query)
    import decimal
    decimalizer = lambda x:(float(x) if type(x) == decimal.Decimal else x)
    return [
        {
            colname:decimalizer(value) for (colname, value) in zip(
                [k[0] for k in column_defs],
                [val for val in row]
            )
        } for row in query.rows()
    ]
