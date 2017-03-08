import os
import csv
from flask import current_app
import subprocess
from .processes import fetch_process, is_running, gen_files_list
from .database import filterfile
from .utils import descriptions, column_filter

def results_get(id):
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
    data = gen_files_list(id, data)
    output = []
    for fileID in process[0]['files']:
        output.append({
            'fileID':fileID,
            'description':process[0]['files'][fileID]['description'],
            'display_name':process[0]['files'][fileID]['display_name'],
            'url':'/api/v1/processes/%d/results/%s'%(id, fileID),
            'size':"%0.3f KB"%(
                os.path.getsize(process[0]['files'][fileID]['fullname'])/1024
            ),
            'rows':int(subprocess.check_output([
                'wc',
                '-l',
                process[0]['files'][fileID]['fullname']
            ]).decode().split()[0])-1,
        })
    return output


def list_input(path = None):
    """Fetches a list of input files from the input directory"""
    data = current_app.config['storage']['loader']()
    if not path:
        path = os.path.join(current_app.config['files']['data-dir'], 'input')
        current_app.config['storage']['manifest'] = []
    output = []
    for entity in sorted(os.listdir(path)):
        fullname = os.path.join(path, entity)
        if os.path.isfile(fullname):
            output.append({
                # 'display_name':entity,
                'name':entity, #not using fullname anymore for security reasons
                'type':'file',
                'fileID':len(current_app.config['storage']['manifest']),
                'description':descriptions(
                    '.'.join(os.path.basename(entity).split('.')[1:])
                ),
            })
            current_app.config['storage']['manifest'].append(fullname)
        elif os.path.isdir(fullname):
            output.append({
                'name':entity, #not using fullname anymore for security reasons
                # 'display_name':entity,
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
                'archive',
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
    data = gen_files_list(id, data)
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
            'size':"%0.3f KB"%(
                os.path.getsize(os.path.join(
                    current_app.config['files']['data-dir'],
                    'archive',
                    entry['display_name']
                ))/1024
            ),
            'rows':int(subprocess.check_output(['wc', '-l', os.path.join(
                current_app.config['files']['data-dir'],
                'archive',
                entry['display_name']
            )]).decode().split()[0])-1,
        } for (key, entry) in data['dropbox'].items()
    ]
