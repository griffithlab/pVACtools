import os
import csv
import re
from flask import current_app
import subprocess
from .processes import fetch_process, is_running
from .database import filterfile
from .utils import descriptions, column_filter, filterdata, sort, fullresponse

# details for each file to be appended to the output of results_get
def resultfile(id, process, fileID):
    return({
        'fileID':int(fileID),
        'description':process[0]['files'][fileID]['description'],
        'display_name':process[0]['files'][fileID]['display_name'],
        'url':'/api/v1/processes/%d/results/%s'%(id, fileID),
        'size':os.path.getsize(process[0]['files'][fileID]['fullname']),
        'rows':int(subprocess.check_output([
            'wc',
            '-l',
            process[0]['files'][fileID]['fullname']
        ]).decode().split()[0])-1,
    })

def results_get(id, type, filters, sorting, page, count):
    """Get the list of result files from a specific pVAC-Seq run"""
    data = current_app.config['storage']['loader']()
    if id == -1:
        return list_dropbox()
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    output = []
    if "default" in type:
         for fileID in process[0]['files']:
            if not (re.search("/tmp/", process[0]['files'][fileID]['display_name'])):
                output.append(resultfile(id,process,fileID))
    elif "all" in type:
        for fileID in process[0]['files']:
            output.append(resultfile(id,process,fileID))
    else:
        for filter in type:
            for fileID in process[0]['files']:
                if (re.search('%s.tsv'%filter, process[0]['files'][fileID]['display_name'])):
                    output.append(resultfile(id,process,fileID))
    return filterdata(output, filters, sorting, page, count)


def list_input(path = None):
    """Fetches a list of input files from the input directory"""
    data = current_app.config['storage']['loader']()
    if not path:
        path = os.path.join(current_app.config['files']['data-dir'], 'input')
        current_app.config['storage']['manifest'] = []
    output = []
    for entity in sorted(os.listdir(path)):
        fullname = os.path.join(path, entity)
        if (fullname[fullname.rfind('/')+1] == '.'):
            continue
        if os.path.isfile(fullname):
            output.append({
                'display_name':entity,
                'type':'file',
                'fileID':len(current_app.config['storage']['manifest']),
                'description':descriptions(
                    '.'.join(os.path.basename(entity).split('.')[1:])
                ),
            })
            current_app.config['storage']['manifest'].append(fullname)
        elif os.path.isdir(fullname):
            output.append({
                'display_name':entity,
                'type':'directory',
                'contents': list_input(fullname)
            })
    return output

def results_getfile(id, fileID, count, page, filters, sort, direction):
    """(DEPRECATED) Read data directly from a specific output file"""
    return filterfile(
        parentID = id,
        fileID = fileID,
        count = count,
        page = page,
        filters = filters,
        sort = sort,
        direction = direction
    )



def results_getcols(id, fileID):
    """Get a mapping of standardized column names -> original column names"""
    data = current_app.config['storage']['loader']()
    if id==-1:
        if str(fileID) not in data['dropbox']:
            return ({
                'code':400,
                'message':'The requested file (%d) does not exist'%fileID,
                'fields':'fileID'
            },400)
        raw_reader = open(
            os.path.join(
                os.path.abspath(current_app.config['files']['data-dir']),
                'dropbox',
                data['dropbox'][str(fileID)]['display_name']
            )
        )
        reader = csv.DictReader(raw_reader, delimiter='\t')
        output = {
            column_filter(field):field for field in reader.fieldnames
        }
        raw_reader.close()
        return output
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code": 400,
                "message": "The requested process (%d) does not exist"%id,
                "fields": "id"
            },400
        )
    if str(fileID) not in process[0]['files']:
        return (
            {
                "code":400,
                "message":"The requested fileID (%d) does not exist for this process (%d)" %(fileID, id),
                "fields":"fileID"
            }, 400
        )
    raw_reader = open(os.path.join(
        process[0]['output'],
        process[0]['files'][str(fileID)]['display_name']
    ))
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = {column_filter(field):field for field in reader.fieldnames}
    raw_reader.close()
    return output

def list_dropbox():
    data = current_app.config['storage']['loader']()
    return [
        {
            'fileID':key,
            'description':entry['description'],
            'display_name':entry['display_name'],
            'url':'/api/v1/processes/-1/results/%s'%(key),
            'size':(
                os.path.getsize(os.path.join(
                    current_app.config['files']['data-dir'],
                    'dropbox',
                    entry['display_name']
                ))
            ),
            'rows':int(subprocess.check_output(['wc', '-l', os.path.join(
                current_app.config['files']['data-dir'],
                'dropbox',
                entry['display_name']
            )]).decode().split()[0])-1,
        } for (key, entry) in data['dropbox'].items()
    ]
